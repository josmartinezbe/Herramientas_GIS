# -*- coding: utf-8 -*-
import arcpy, os, csv, sys

arcpy.env.overwriteOutput = True

class Toolbox(object):
    def __init__(self):
        self.label = "Priorización de Predios (Vectorial, WLC)"
        self.alias = "priorizacion_predios"
        self.tools = [PriorizarPredios]

class PriorizarPredios(object):
    def __init__(self):
        self.label = "Calcular Índice WLC por Predio"
        self.description = (
            "Proyecta insumos a un SR objetivo (CTM12 recomendado) y calcula, por predio, "
            "Pct_<criterio>, Score_<criterio>, Indice_WLC, Indice_Efic y Rank_WLC "
            "a partir de criterios vectoriales con un campo 'Valor' y un CSV de pesos."
        )
        self.canRunInBackground = True

    # --------------------------
    # Parámetros
    # --------------------------
    def getParameterInfo(self):
        p_gdb = arcpy.Parameter("Geodatabase (GDB)", "gdb",
                                "DEWorkspace", "Required", "Input")

        p_fds_criterios = arcpy.Parameter("Feature Dataset de Criterios",
                                "fds_criterios", "DEFeatureDataset",
                                "Required", "Input")
        p_fds_criterios.parameterDependencies = [p_gdb.name]

        p_fds_predios = arcpy.Parameter("Feature Dataset de Predios",
                                "fds_predios", "DEFeatureDataset",
                                "Required", "Input")
        p_fds_predios.parameterDependencies = [p_gdb.name]

        p_fc_predios = arcpy.Parameter("Feature Class de Predios",
                                "fc_predios", "DEFeatureClass",
                                "Required", "Input")
        p_fc_predios.parameterDependencies = [p_fds_predios.name]

        p_campo_id = arcpy.Parameter("Campo ID de Predio",
                                "campo_id_predio", "Field",
                                "Required", "Input")
        p_campo_id.parameterDependencies = [p_fc_predios.name]

        p_campo_valor = arcpy.Parameter("Nombre del campo de Valor en criterios",
                                "campo_valor", "GPString",
                                "Required", "Input")
        p_campo_valor.value = "Valor"

        p_csv_pesos = arcpy.Parameter(
            "CSV de Pesos y valor_max (criterio,peso,valor_max)",
            "csv_pesos", "DEFile", "Required", "Input")
        p_csv_pesos.filter.list = ["csv"]

        p_fds_salida = arcpy.Parameter("Feature Dataset de Salida (se creará si no existe)",
                                "fds_salida", "DEFeatureDataset",
                                "Required", "Input")
        p_fds_salida.parameterDependencies = [p_gdb.name]

        p_nombre_salida = arcpy.Parameter("Nombre del Feature Class de Salida",
                                "fc_salida", "GPString",
                                "Required", "Input")
        p_nombre_salida.value = "Predios_Prioridad"

        p_calc_efic = arcpy.Parameter("Calcular Índice de Eficiencia (Indice_WLC / Area_ha)",
                                "calc_efic", "GPBoolean",
                                "Optional", "Input")
        p_calc_efic.value = True

        p_umbral_area = arcpy.Parameter("Umbral de Área (ha) para Ranking",
                                "umbral_area", "GPDouble",
                                "Optional", "Input")
        p_umbral_area.value = 100.0

        # SR de salida opcional: si no se elige, intentará 9377 y luego 3116
        p_sr_out = arcpy.Parameter(
            displayName="Sistema de referencia de salida (opcional; recom.: CTM12). Puedes seleccionar un SR o dar un .prj",
            name="sr_out",
            datatype="GPCoordinateSystem",
            parameterType="Optional",
            direction="Input"
        )

        return [p_gdb, p_fds_criterios, p_fds_predios, p_fc_predios, p_campo_id,
                p_campo_valor, p_csv_pesos, p_fds_salida, p_nombre_salida,
                p_calc_efic, p_umbral_area, p_sr_out]

    def updateMessages(self, params):
        gdb, fds_crit, fds_pred, fc_pred, _, _, csv_pesos, fds_out = \
            [p.valueAsText for p in params[:8]]
        if gdb and not arcpy.Exists(gdb):
            params[0].setErrorMessage("La GDB no existe.")
        if fds_crit and not arcpy.Exists(fds_crit):
            params[1].setErrorMessage("El Feature Dataset de criterios no existe.")
        if fds_pred and not arcpy.Exists(fds_pred):
            params[2].setErrorMessage("El Feature Dataset de predios no existe.")
        if fc_pred and not arcpy.Exists(fc_pred):
            params[3].setErrorMessage("La Feature Class de predios no existe.")
        if csv_pesos and not os.path.isfile(csv_pesos):
            params[6].setErrorMessage("El CSV de pesos no existe.")

    # --------------------------
    # Utilidades
    # --------------------------
    def _get_sr_out(self, sr_param_text):
        """
        Devuelve un SpatialReference válido.
        Prioridad: parámetro de usuario -> 9377 -> 3116 -> error instructivo.
        """
        if sr_param_text:
            # puede ser nombre, WKID o ruta .prj
            try:
                return arcpy.SpatialReference(sr_param_text)
            except:
                pass
            if os.path.isfile(sr_param_text) and sr_param_text.lower().endswith(".prj"):
                return arcpy.SpatialReference(sr_param_text)
            raise RuntimeError(u"No se pudo cargar el SR indicado: {}".format(sr_param_text))

        # 9377 (CTM12)
        try:
            sr = arcpy.SpatialReference(9377)
            if sr and sr.factoryCode not in (0, None):
                return sr
        except:
            pass

        # 3116 (MAGNA-SIRGAS / Colombia Bogotá) — común en 10.x
        try:
            sr = arcpy.SpatialReference(3116)
            if sr and sr.factoryCode not in (0, None):
                return sr
        except:
            pass

        raise RuntimeError(
            u"No se encontró un SR válido. Selecciona uno en el parámetro "
            u"“Sistema de referencia de salida” o proporciona un .prj."
        )

    def _add_field_if_absent(self, fc, name, ftype="DOUBLE"):
        if name not in [f.name for f in arcpy.ListFields(fc)]:
            arcpy.AddField_management(fc, name, ftype)

    def _compute_area_ha(self, fc):
        self._add_field_if_absent(fc, "Area_ha", "DOUBLE")
        with arcpy.da.UpdateCursor(fc, ["SHAPE@AREA", "Area_ha"]) as cur:
            for a_m2, _ in cur:
                cur.updateRow([a_m2, (a_m2 or 0.0)/10000.0])

    def _ensure_fds(self, gdb_path, fds_name, sr):
        fds_path = os.path.join(gdb_path, fds_name)
        if not arcpy.Exists(fds_path):
            arcpy.CreateFeatureDataset_management(gdb_path, fds_name, sr)
        return fds_path

    def _project_to_sr(self, in_fc, out_fc, target_sr):
        sr_in = arcpy.Describe(in_fc).spatialReference
        if (not sr_in) or sr_in.name in (u"Unknown", u"Undefined") or sr_in.factoryCode in (0, None):
            raise RuntimeError(
                u"'{}' no tiene SR definido. Usa 'Define Projection' al SR verdadero antes de ejecutar."
                .format(in_fc)
            )
        if sr_in.factoryCode == target_sr.factoryCode and sr_in.name == target_sr.name:
            arcpy.CopyFeatures_management(in_fc, out_fc)
        else:
            arcpy.Project_management(in_fc, out_fc, target_sr)
        return out_fc

    def _max_valor_in_fc(self, fc, campo_valor):
        vmax = 0.0
        with arcpy.da.SearchCursor(fc, [campo_valor]) as cur:
            for (v,) in cur:
                try:
                    vmax = max(vmax, float(v) if v is not None else 0.0)
                except:
                    pass
        return vmax if vmax > 0 else 1.0

    def _read_pesos(self, csv_path):
        """Lectura robusta (Py2/Py3) de CSV: criterio,peso,valor_max"""
        def to_text(b):
            if b is None: return u""
            if sys.version_info[0] < 3:
                if isinstance(b, unicode):  # noqa
                    return b
                try:
                    return b.decode('utf-8')
                except:
                    try:
                        return b.decode('latin-1')
                    except:
                        return unicode(str(b), errors='ignore')  # noqa
            else:
                return str(b)

        def to_float(x, default=0.0):
            s = to_text(x).strip().replace(',', '.')
            try:
                return float(s) if s != u'' else default
            except:
                return default

        pesos, vmaxs = {}, {}
        if sys.version_info[0] < 3:
            f = open(csv_path, 'rb')
            try:
                rdr = csv.DictReader(f)
                if rdr.fieldnames and rdr.fieldnames[0].startswith('\xef\xbb\xbf'):
                    rdr.fieldnames[0] = rdr.fieldnames[0][3:]
                for row in rdr:
                    crit = to_text(row.get('criterio')).strip() or to_text(row.get('\xef\xbb\xbfcriterio')).strip()
                    if not crit: continue
                    pesos[crit] = to_float(row.get('peso'), 0.0)
                    vmax_raw = row.get('valor_max')
                    vmaxs[crit] = to_float(vmax_raw, None) if (vmax_raw is not None and to_text(vmax_raw).strip() != u'') else None
            finally:
                f.close()
        else:
            with open(csv_path, 'r', encoding='utf-8-sig', newline='') as f:
                rdr = csv.DictReader(f)
                for row in rdr:
                    crit = (row.get('criterio') or '').strip()
                    if not crit: continue
                    pesos[crit] = to_float(row.get('peso'), 0.0)
                    vmax_txt = row.get('valor_max')
                    vmaxs[crit] = to_float(vmax_txt, None) if (vmax_txt is not None and vmax_txt.strip() != '') else None

        return pesos, vmaxs

    def _intersect_score(self, predios_fc, criterio_fc, id_field, campo_valor, valor_max):
        ix = arcpy.CreateUniqueName("ix_" + os.path.basename(criterio_fc), "in_memory")
        arcpy.Intersect_analysis([predios_fc, criterio_fc], ix, "ALL", "", "INPUT")

        self._add_field_if_absent(ix, "Area_ix_ha", "DOUBLE")
        with arcpy.da.UpdateCursor(ix, ["SHAPE@AREA", "Area_ix_ha"]) as cur:
            for a_m2, _ in cur:
                cur.updateRow([a_m2, (a_m2 or 0.0)/10000.0])

        area_pred = {}
        with arcpy.da.SearchCursor(predios_fc, [id_field, "Area_ha"]) as cur:
            for pid, a in cur:
                area_pred[pid] = a or 0.0

        scores, covers = {}, {}
        fields = [id_field, "Area_ix_ha", campo_valor]
        with arcpy.da.SearchCursor(ix, fields) as cur:
            for pid, a_ix, val in cur:
                if pid is None or a_ix is None: continue
                a_p = area_pred.get(pid, 0.0)
                if a_p <= 0: continue
                pct = a_ix / a_p
                vnorm = (float(val)/float(valor_max)) if valor_max else 0.0
                scores[pid] = scores.get(pid, 0.0) + (pct * vnorm)
                covers[pid] = covers.get(pid, 0.0) + pct

        arcpy.Delete_management(ix)
        return scores, covers

    # --------------------------
    # Ejecutar
    # --------------------------
    def execute(self, params, messages):
        # Parámetros
        GDB            = params[0].valueAsText
        FDS_CRITERIOS  = params[1].valueAsText
        FDS_PREDIOS    = params[2].valueAsText
        FC_PREDIOS     = params[3].valueAsText
        CAMPO_ID_PRED  = params[4].valueAsText
        CAMPO_VALOR    = params[5].valueAsText
        CSV_PESOS      = params[6].valueAsText
        FDS_SALIDA     = params[7].valueAsText
        FC_SALIDA_NOM  = params[8].valueAsText
        CALC_EFIC      = bool(params[9].value)
        UMBRAL_AREA_HA = float(params[10].value)
        SR_PARAM_TXT   = params[11].valueAsText  # SR opcional del usuario

        arcpy.AddMessage("Leyendo pesos...")
        PESOS, VMAXS = self._read_pesos(CSV_PESOS)
        arcpy.AddMessage(u"Criterios en CSV: {}".format(", ".join(sorted(PESOS.keys()))))

        # 1) SR de salida (con fallback 9377 -> 3116 si no se define)
        SR_OUT = self._get_sr_out(SR_PARAM_TXT)
        arcpy.AddMessage(u"SR de salida: {} (WKID: {})".format(SR_OUT.name, SR_OUT.factoryCode))

        # 2) FDS temporal en SR_OUT
        fds_tmp_name = "_TMP_SR"
        fds_tmp = self._ensure_fds(GDB, fds_tmp_name, SR_OUT)

        # 3) Predios -> SR_OUT
        pred_out = os.path.join(fds_tmp, "Predios_SR")
        self._project_to_sr(FC_PREDIOS, pred_out, SR_OUT)
        self._compute_area_ha(pred_out)

        # 4) Criterios -> SR_OUT
        arcpy.env.workspace = FDS_CRITERIOS
        criterios_src = arcpy.ListFeatureClasses()
        if not criterios_src:
            raise RuntimeError("No se encontraron Feature Classes de criterios en {}".format(FDS_CRITERIOS))

        criterios_proj = []
        for fc_name in criterios_src:
            src = os.path.join(FDS_CRITERIOS, fc_name)
            dst = os.path.join(fds_tmp, "{}_SR".format(fc_name))
            self._project_to_sr(src, dst, SR_OUT)
            criterios_proj.append((fc_name, dst))

        # 5) Preparar salida
        fds_out_dir = os.path.dirname(FDS_SALIDA)
        fds_out_name = os.path.basename(FDS_SALIDA)
        if not arcpy.Exists(FDS_SALIDA):
            arcpy.CreateFeatureDataset_management(fds_out_dir, fds_out_name, SR_OUT)
        out_fc = os.path.join(FDS_SALIDA, FC_SALIDA_NOM)
        if arcpy.Exists(out_fc):
            arcpy.Delete_management(out_fc)
        arcpy.CopyFeatures_management(pred_out, out_fc)

        # 6) Intersectar y calcular por criterio
        score_fields = []
        for fc_name, fc_proj in criterios_proj:
            crit_key = fc_name  # debe existir en el CSV (columna 'criterio')
            if crit_key not in PESOS:
                arcpy.AddWarning(u"Omitido: '{}' (sin peso en CSV)".format(crit_key))
                continue

            peso = float(PESOS[crit_key])
            vmax = VMAXS.get(crit_key)
            if vmax is None:
                vmax = self._max_valor_in_fc(fc_proj, CAMPO_VALOR)

            scores, covers = self._intersect_score(pred_out, fc_proj, CAMPO_ID_PRED, CAMPO_VALOR, vmax)

            # nombres de campo validados (para evitar truncamientos/caracteres inválidos)
            raw_score = u"Score_{}".format(crit_key)
            raw_pct   = u"Pct_{}".format(crit_key)
            f_score = arcpy.ValidateFieldName(raw_score, out_fc)
            f_pct   = arcpy.ValidateFieldName(raw_pct, out_fc)

            self._add_field_if_absent(out_fc, f_score, "DOUBLE")
            self._add_field_if_absent(out_fc, f_pct,   "DOUBLE")

            with arcpy.da.UpdateCursor(out_fc, [CAMPO_ID_PRED, f_score, f_pct]) as cur:
                for pid, _, _ in cur:
                    cur.updateRow([pid, scores.get(pid, 0.0), covers.get(pid, 0.0)])

            score_fields.append((f_score, peso))
            arcpy.AddMessage(u"Incluido criterio: '{}'  → campo: {}  peso: {}".format(crit_key, f_score, peso))

        arcpy.AddMessage(u"Criterios realmente incluidos: {}".format(", ".join([f for (f, _) in score_fields])))

        # 7) Índice WLC (robusto)
        if not score_fields:
            raise RuntimeError(
                u"No se agregó ningún criterio desde el CSV. "
                u"Revisa que la columna 'criterio' coincida EXACTAMENTE con los nombres de los FC "
                u"en el FDS de criterios."
            )

        self._add_field_if_absent(out_fc, "Indice_WLC", "DOUBLE")
        fields_wlc = [fld for (fld, _) in score_fields] + ["Indice_WLC"]

        with arcpy.da.UpdateCursor(out_fc, fields_wlc) as cur:
            idx2w = {i: w for i, (_, w) in enumerate(score_fields)}
            for row in cur:
                total = 0.0
                for i, w in idx2w.items():
                    try:
                        v = row[i]
                    except IndexError:
                        v = 0.0
                    total += (0.0 if v is None else float(v)) * w
                row[-1] = total
                cur.updateRow(row)

        # 8) Eficiencia y ranking
        if bool(params[9].value):
            self._add_field_if_absent(out_fc, "Indice_Efic", "DOUBLE")
            with arcpy.da.UpdateCursor(out_fc, ["Area_ha", "Indice_WLC", "Indice_Efic"]) as cur:
                for a_ha, idx, _ in cur:
                    eff = (idx / a_ha) if (a_ha and a_ha > 0) else None
                    cur.updateRow([a_ha, idx, eff])

        self._add_field_if_absent(out_fc, "Rank_WLC", "LONG")
        rows = []
        with arcpy.da.SearchCursor(out_fc, [CAMPO_ID_PRED, "Area_ha", "Indice_WLC"]) as cur:
            for pid, a_ha, idx in cur:
                if a_ha and a_ha >= UMBRAL_AREA_HA:
                    rows.append((pid, a_ha, idx))
        rows.sort(key=lambda r: (r[2] if r[2] is not None else -1), reverse=True)
        rank_map = {pid: i+1 for i,(pid,_,_) in enumerate(rows)}

        with arcpy.da.UpdateCursor(out_fc, [CAMPO_ID_PRED, "Area_ha", "Rank_WLC"]) as cur:
            for pid, a_ha, _ in cur:
                cur.updateRow([pid, a_ha, rank_map.get(pid, None)])

        arcpy.AddMessage(u"✅ Listo. Salida: {}".format(out_fc))
        arcpy.AddMessage(u"ℹ️ Insumos proyectados temporalmente a: {}".format(os.path.join(GDB, "_TMP_SR")))

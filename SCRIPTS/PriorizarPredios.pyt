# -*- coding: utf-8 -*-
# PriorizarPredios.pyt - ArcMap 10.8 (Python 2.7)
import arcpy, os, sys, re, unicodedata
arcpy.env.overwriteOutput = True

# ---------------- utilidades generales ----------------
def to_text(b):
    if b is None:
        return u""
    if sys.version_info[0] < 3:
        try:
            return b if isinstance(b, unicode) else b.decode('utf-8')
        except:
            try:
                return b.decode('latin-1')
            except:
                return unicode(str(b), errors='ignore')
    return unicode(b)

def strip_accents(s):
    s = to_text(s)
    return u"".join([c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn'])

def norm_name(s):
    """Normaliza nombres para emparejar FC <-> criterio (sin tildes, espacios, sufijos)."""
    if not s:
        return u""
    s = strip_accents(to_text(s)).strip().upper()
    for suf in ["_AOI_WGS84", "_WGS84", "_AOI", "_SR"]:
        if s.endswith(suf):
            s = s[:-len(suf)]
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s

def gdb_of_path(path_in_gdb):
    if not path_in_gdb:
        return None
    p = path_in_gdb.replace("\\", "/")
    u = p.upper()
    i = u.rfind(".GDB/")
    if i == -1:
        i = u.rfind(".GDB")
    return p[:i+4] if i != -1 else os.path.dirname(p)

def list_fcs_in_fds(fds_path):
    fcs = []
    if not fds_path or not arcpy.Exists(fds_path):
        return fcs
    for dp, dns, fns in arcpy.da.Walk(fds_path, datatype="FeatureClass"):
        for fn in fns:
            try:
                nm = arcpy.Describe(os.path.join(dp, fn)).baseName
            except:
                nm = os.path.basename(fn)
            if nm.upper().startswith("PREDIOS_"):
                continue
            fcs.append(nm)
    return sorted(set(fcs))

# ----------------- toolbox -----------------
class Toolbox(object):
    def __init__(self):
        self.label = "Priorizacion de Predios (WLC con costos que restan)"
        self.alias = "priorizacion_predios"
        self.tools = [PriorizarPredios]

class PriorizarPredios(object):
    def __init__(self):
        self.label = "Calcular Indice WLC (beneficio suma, costo resta)"
        self.description = ("Lee parametros desde tbl_config. Intersecta, repara geometria, "
                            "convierte a singlepart y disuelve por predio. Calcula Pct_, Score_, Peso_, "
                            "Vmin_, Vmax_, Comp_, Benef_Total, Cost_Total e Indice_Neto.")
        self.canRunInBackground = True

    # ---------- parametros ----------
    def getParameterInfo(self):
        p_gdb = arcpy.Parameter(
            displayName="Geodatabase (GDB)",
            name="gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input"
        )

        p_fds_criterios = arcpy.Parameter(
            displayName="Feature Dataset de Criterios",
            name="fds_criterios",
            datatype="DEFeatureDataset",
            parameterType="Required",
            direction="Input"
        )
        p_fds_criterios.parameterDependencies = [p_gdb.name]

        p_fds_predios = arcpy.Parameter(
            displayName="Feature Dataset de Predios",
            name="fds_predios",
            datatype="DEFeatureDataset",
            parameterType="Required",
            direction="Input"
        )
        p_fds_predios.parameterDependencies = [p_gdb.name]

        p_fc_predios = arcpy.Parameter(
            displayName="Feature Class de Predios",
            name="fc_predios",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )
        p_fc_predios.parameterDependencies = [p_fds_predios.name]

        p_campo_id = arcpy.Parameter(
            displayName="Campo ID de Predio",
            name="campo_id_predio",
            datatype="Field",
            parameterType="Required",
            direction="Input"
        )
        p_campo_id.parameterDependencies = [p_fc_predios.name]

        p_campo_valor = arcpy.Parameter(
            displayName="Campo 'Valor' en criterios",
            name="campo_valor",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        p_campo_valor.value = "Valor"

        p_tbl_config = arcpy.Parameter(
            displayName="Tabla de configuracion (tbl_config)",
            name="tbl_config",
            datatype="DETable",
            parameterType="Required",
            direction="Input"
        )

        # En ArcMap/10.8 usar GPSpatialReference
        p_sr_out = arcpy.Parameter(
            displayName="SR de salida - 3116 recomendado",
            name="sr_out",
            datatype="GPSpatialReference",
            parameterType="opcional",
            direction="Input"
        )

        p_fds_salida = arcpy.Parameter(
            displayName="Feature Dataset de Salida",
            name="fds_salida",
            datatype="DEFeatureDataset",
            parameterType="Required",
            direction="Input"
        )
        p_fds_salida.parameterDependencies = [p_gdb.name]

        p_nombre_salida = arcpy.Parameter(
            displayName="Nombre FC de salida",
            name="fc_salida",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        p_nombre_salida.value = "Predios_Prioridad"

        return [p_gdb, p_fds_criterios, p_fds_predios, p_fc_predios, p_campo_id,
                p_campo_valor, p_tbl_config, p_sr_out, p_fds_salida, p_nombre_salida]

    def updateMessages(self, params):
        idxs_exist = [0, 1, 2, 3, 6, 8]
        for p in idxs_exist:
            if params[p].valueAsText and not arcpy.Exists(params[p].valueAsText):
                params[p].setErrorMessage("La ruta no existe o no es accesible.")
        if params[8].valueAsText and not arcpy.Exists(params[8].valueAsText):
            params[8].setWarningMessage("El FDS de salida se creara en la ejecucion.")

    # ---------- helpers espaciales ----------
    def get_sr_out(self, sr_param_obj, fc_predios, fds_criterios):
        try:
            if sr_param_obj:
                fc = getattr(sr_param_obj, "factoryCode", None)
                nm = getattr(sr_param_obj, "name", None)
                if fc not in (None, 0):
                    return arcpy.SpatialReference(fc)
                if nm and nm not in ("Unknown", "Undefined"):
                    return sr_param_obj
        except:
            pass
        try:
            sr_p = arcpy.Describe(fc_predios).spatialReference
            if sr_p and sr_p.name not in ("Unknown","Undefined") and sr_p.factoryCode not in (0, None):
                return sr_p
        except:
            pass
        try:
            for dp, dns, fns in arcpy.da.Walk(fds_criterios, datatype="FeatureClass"):
                for fn in fns:
                    sr_c = arcpy.Describe(os.path.join(dp, fn)).spatialReference
                    if sr_c and sr_c.name not in ("Unknown","Undefined") and sr_c.factoryCode not in (0, None):
                        return sr_c
        except:
            pass
        return arcpy.SpatialReference(3116)

    def ensure_fds(self, gdb_path, fds_name, sr):
        """Crea (si hace falta) y devuelve la ruta del Feature Dataset dentro de la GDB."""
        if fds_name.upper().endswith(".GDB") or ".GDB" in fds_name.upper():
            gdb_path = (gdb_of_path(fds_name) or gdb_path)
            fds_name = os.path.basename(fds_name)
        fds_path = os.path.join(gdb_path, fds_name)
        if not arcpy.Exists(fds_path):
            arcpy.CreateFeatureDataset_management(gdb_path, fds_name, sr)
        else:
            try:
                sr_exist = arcpy.Describe(fds_path).spatialReference
                if (sr_exist and getattr(sr_exist, "factoryCode", None) not in (None, 0) and
                    getattr(sr, "factoryCode", None) not in (None, 0) and
                    sr_exist.factoryCode != sr.factoryCode):
                    raise RuntimeError(
                        "El FDS '{}' ya existe con SR {} ({}), distinto al SR destino {} ({})."
                        .format(fds_name, sr_exist.name, sr_exist.factoryCode, sr.name, sr.factoryCode)
                    )
            except:
                pass
        return fds_path

    def add_field_if_absent(self, fc, name, ftype="DOUBLE"):
        fields = [f.name for f in arcpy.ListFields(fc)]
        if name not in fields:
            arcpy.AddField_management(fc, name, ftype)

    def compute_area_ha(self, fc):
        self.add_field_if_absent(fc, "Area_ha", "DOUBLE")
        with arcpy.da.UpdateCursor(fc, ["SHAPE@AREA", "Area_ha"]) as cur:
            for a_m2, _ in cur:
                cur.updateRow([a_m2, (a_m2 or 0.0) / 10000.0])

    def pick_geo_transform(sr_in, sr_out):
        try:
            cands = arcpy.ListTransformations(sr_in, sr_out)
            if not cands:
                return None
            # prioriza transformaciones MAGNA/Colombia si existen
            pref = [t for t in cands if ("MAGNA" in t.upper() or "COLOMB" in t.upper())]
            return pref[0] if pref else cands[0]
        except:
            return None
    def project_to_sr(self, in_fc, out_fc, target_sr):
        sr_in = arcpy.Describe(in_fc).spatialReference
        # solo invalidar si es Unknown/Undefined
        if (not sr_in) or sr_in.name in ("Unknown", "Undefined"):
            raise RuntimeError("La capa '{}' no tiene SR definido. Use 'Define Projection' y reintente."
                            .format(in_fc))
        # si ya está en el SR destino, copiar
        if (sr_in.factoryCode == target_sr.factoryCode and
            sr_in.name == target_sr.name):
            arcpy.CopyFeatures_management(in_fc, out_fc)
            return out_fc

        # elegir transformación si hace falta (p.ej., GCS distintos)
        try:
            g_in  = getattr(sr_in,  "GCS", None) or getattr(sr_in,  "geographicCoordinateSystem", None)
            g_out = getattr(target_sr, "GCS", None) or getattr(target_sr, "geographicCoordinateSystem", None)
            need_gt = (g_in and g_out and g_in.name != g_out.name)
        except:
            need_gt = False

        gt = pick_geo_transform(sr_in, target_sr) if need_gt else None
        if gt:
            arcpy.Project_management(in_fc, out_fc, target_sr, gt)
        else:
            arcpy.Project_management(in_fc, out_fc, target_sr)
        return out_fc
    # ---------- lectura de configuracion ----------
        # ---------- lectura de configuracion (robusta a columnas faltantes) ----------
    def read_config(self, tbl_path):
        """
        Lee tbl_config aceptando columnas opcionales ausentes.
        Columnas esperadas (nombres flexibles):
          - Criterio / Criterio_norm / FC / Nombre     (obligatorio)
          - Peso / Weight                              (obligatorio)
          - Metodo                                     (obligatorio)
          - Polarity (opcional; default 'benefit')
          - Vmin, Vmax (opcionales; default None)
          - ConstVal / Const (opcional; default None)
        Normaliza los pesos globalmente (suman 1) sobre criterios con peso > 0.
        """
        campos = [f.name for f in arcpy.ListFields(tbl_path)]
        cmap = {c.upper(): c for c in campos}

        def pick(*cands):
            for c in cands:
                if c.upper() in cmap:
                    return cmap[c.upper()]
            return None

        c_nom  = pick("Criterio", "Criterio_norm", "FC", "Nombre")
        c_w    = pick("Peso", "Weight")
        c_met  = pick("Metodo", "Method")
        c_pol  = pick("Polarity")
        c_vmin = pick("Vmin", "Min")
        c_vmax = pick("Vmax", "Max")
        c_cst  = pick("ConstVal", "Const")

        # Validar obligatorios
        need = [c_nom, c_w, c_met]
        if any(v is None for v in need):
            raise RuntimeError(
                "tbl_config debe tener al menos: Criterio, Peso y Metodo."
            )

        # Construir la lista de campos REALMENTE presentes
        sel_fields = [c_nom, c_w, c_met]
        for c in (c_pol, c_vmin, c_vmax, c_cst):
            if c: sel_fields.append(c)

        idx = {}  # mapa nombre->indice en el cursor
        for i, fn in enumerate(sel_fields):
            idx[fn] = i

        cfg = {}
        with arcpy.da.SearchCursor(tbl_path, sel_fields) as cur:
            for row in cur:
                nom = row[idx[c_nom]]
                k = norm_name(nom)
                if not k:
                    continue

                # Peso
                try:
                    w = float(row[idx[c_w]])
                except:
                    w = 0.0

                # Metodo (minmax/max/binario/const)
                met = to_text(row[idx[c_met]]).lower() if row[idx[c_met]] not in (None, "") else "minmax"

                # Opcionales con defaults
                pol  = "benefit"
                vmin = None
                vmax = None
                cst  = None

                if c_pol in idx:
                    pol = to_text(row[idx[c_pol]]).lower() if row[idx[c_pol]] not in (None, "") else "benefit"
                if c_vmin in idx:
                    try: vmin = None if row[idx[c_vmin]] in (None, "") else float(row[idx[c_vmin]])
                    except: vmin = None
                if c_vmax in idx:
                    try: vmax = None if row[idx[c_vmax]] in (None, "") else float(row[idx[c_vmax]])
                    except: vmax = None
                if c_cst in idx:
                    try: cst = None if row[idx[c_cst]]  in (None, "") else float(row[idx[c_cst]])
                    except: cst = None

                if w > 0:
                    cfg[k] = {
                        "peso": w,
                        "metodo": met,
                        "polarity": pol,
                        "vmin": vmin,
                        "vmax": vmax,
                        "const": cst,
                        "nombre": nom
                    }

        # Normalizar pesos (suman 1)
        s = sum(v["peso"] for v in cfg.values())
        if s > 0:
            for k2 in cfg:
                cfg[k2]["peso"] = cfg[k2]["peso"] / s

        return cfg


    # ---------- interseccion + singlepart + dissolve ----------
    def intersect_score(self, predios_fc, criterio_fc, id_field, campo_valor, cfg_i):
        """
        Intersecta predios x criterio, repara, singlepart, disuelve por predio
        y devuelve dicts: scores (0..1) / covers (0..1) por ID de predio.

        NOTA: Aqui NO se invierte por 'polarity'; el score siempre representa
        "intensidad" (0..1). La penalizacion de 'cost' se aplica como signo
        al componenete ponderado (Comp_i) en el paso de escritura.
        """
        # 0) copiar criterio y crear VALOR_CRIT (DOUBLE)
        tmp_crit = arcpy.CreateUniqueName("crit_tmp", "in_memory")
        arcpy.CopyFeatures_management(criterio_fc, tmp_crit)
        nombres_tmp = [f.name.upper() for f in arcpy.ListFields(tmp_crit)]
        if "VALOR_CRIT" not in nombres_tmp:
            arcpy.AddField_management(tmp_crit, "VALOR_CRIT", "DOUBLE")

        # campo fuente
        campo_src = None
        for f in arcpy.ListFields(tmp_crit):
            if f.name.upper() == campo_valor.upper():
                campo_src = f
                break
        if not campo_src:
            arcpy.Delete_management(tmp_crit)
            raise RuntimeError("El campo '{}' no existe en {}".format(campo_valor, criterio_fc))

        # rellenar VALOR_CRIT
        with arcpy.da.UpdateCursor(tmp_crit, ["VALOR_CRIT", campo_src.name]) as cur:
            for vcrit, vsrc in cur:
                try:
                    v = float(vsrc) if vsrc not in (None, "") else None
                except:
                    v = None
                cur.updateRow([v, vsrc])

        # 1) intersect y reparar
        tmp_ix = arcpy.CreateUniqueName("ix_", "in_memory")
        arcpy.Intersect_analysis([predios_fc, tmp_crit], tmp_ix, "ALL", "", "INPUT")
        try:
            arcpy.RepairGeometry_management(tmp_ix, "DELETE_NULL")
        except:
            pass

        # 2) multipart -> singlepart
        tmp_sp = arcpy.CreateUniqueName("ix_sp", "in_memory")
        arcpy.MultipartToSinglepart_management(tmp_ix, tmp_sp)

        # 3) metricas por segmento
        self.add_field_if_absent(tmp_sp, "Area_ix_ha", "DOUBLE")
        self.add_field_if_absent(tmp_sp, "A_vnorm", "DOUBLE")

        met   = cfg_i.get("metodo","minmax")
        vmin  = cfg_i.get("vmin", 0.0)
        vmax  = cfg_i.get("vmax", 1.0)
        constv= cfg_i.get("const", None)

        # asegurar vmin/vmax segun metodo
        if met in ("minmax","max"):
            # calcular vmax si falta
            if vmax in (None, "", 0):
                vmax = 0.0
                with arcpy.da.SearchCursor(tmp_crit, ["VALOR_CRIT"]) as c:
                    for (v,) in c:
                        try:
                            vmax = max(vmax, float(v) if v is not None else 0.0)
                        except:
                            pass
                if vmax == 0:
                    vmax = 1.0
            # calcular vmin si es minmax y falta o degenera
            if met == "minmax" and (vmin is None or vmin == vmax):
                vmin2 = None
                with arcpy.da.SearchCursor(tmp_crit, ["VALOR_CRIT"]) as c:
                    for (v,) in c:
                        try:
                            fv = float(v)
                            vmin2 = fv if vmin2 is None else min(vmin2, fv)
                        except:
                            pass
                if vmin2 is None or vmin2 == vmax:
                    vmin = 0.0
                else:
                    vmin = vmin2

        den = (vmax - vmin) if met == "minmax" else (vmax if vmax else 1.0)

        with arcpy.da.UpdateCursor(tmp_sp, ["SHAPE@AREA","Area_ix_ha","VALOR_CRIT","A_vnorm"]) as cur:
            for a_m2, _, val, _ in cur:
                area_ha = (a_m2 or 0.0) / 10000.0

                # normalizacion por segmento (sin invertir por 'cost')
                if met == "const":
                    vnorm = float(constv or 0.0)

                elif met == "binario":
                    vnorm = 1.0 if (val not in (None, "") and float(val) > 0) else 0.0

                elif met == "max":
                    vnorm = 0.0 if den <= 0 else (0.0 if val in (None, "") else float(val) / den)

                else:  # minmax (por defecto)
                    x = 0.0 if val in (None, "") else float(val)
                    vnorm = 0.0 if den <= 0 else (x - vmin) / den

                # clamp 0..1 y acumular area*vnorm
                if vnorm < 0.0: vnorm = 0.0
                if vnorm > 1.0: vnorm = 1.0

                cur.updateRow([a_m2, area_ha, val, area_ha * vnorm])

        # 4) disolver por predio
        tmp_ds = arcpy.CreateUniqueName("ix_ds", "in_memory")
        arcpy.Dissolve_management(tmp_sp, tmp_ds, [id_field],
                                  [["Area_ix_ha","SUM"], ["A_vnorm","SUM"]])

        # 5) areas de predio y metricas finales (cap score <= pct <= 1)
        area_pred = {}
        with arcpy.da.SearchCursor(predios_fc, [id_field, "Area_ha"]) as cur:
            for pid, a in cur:
                area_pred[pid] = a or 0.0

        scores, covers = {}, {}
        with arcpy.da.SearchCursor(tmp_ds, [id_field, "SUM_Area_ix_ha", "SUM_A_vnorm"]) as cur:
            for pid, a_ix_ha, a_vn in cur:
                a_p = area_pred.get(pid, 0.0)
                if a_p <= 0:
                    continue
                pct = (a_ix_ha or 0.0) / a_p
                if pct < 0.0: pct = 0.0
                if pct > 1.0: pct = 1.0

                sc  = (a_vn or 0.0) / a_p
                if sc < 0.0: sc = 0.0
                if sc > pct: sc = pct
                if sc > 1.0: sc = 1.0

                scores[pid] = sc
                covers[pid] = pct

        # limpieza
        for t in [tmp_ix, tmp_sp, tmp_ds, tmp_crit]:
            try:
                arcpy.Delete_management(t)
            except:
                pass

        return scores, covers

    # ---------- ejecutar ----------
    def execute(self, params, messages):
        (GDB, FDS_CRITERIOS, FDS_PREDIOS, FC_PREDIOS, CAMPO_ID_PRED,
         CAMPO_VALOR, TBL_CONFIG, SR_PARAM, FDS_SALIDA, FC_SALIDA_NOM) = [p.valueAsText for p in params]

        SR_OUT = arcpy.SpatialReference(3116)  # MAGNA-SIRGAS / Colombia Bogota (EPSG:3116)
        arcpy.env.outputCoordinateSystem = SR_OUT

        # 1) leer configuracion
        cfg = self.read_config(TBL_CONFIG)
        if not cfg:
            raise RuntimeError("No se encontraron criterios con peso > 0 en tbl_config.")
        arcpy.AddMessage("Criterios cargados: {}".format(", ".join(sorted([to_text(v["nombre"]) for v in cfg.values()]))))

        # 2) preparar FDS temporal y predios reproyectados
        out_gdb = gdb_of_path(FDS_SALIDA) or GDB
        fds_tmp = self.ensure_fds(out_gdb, "_TMP_SR", SR_OUT)

        pred_out = os.path.join(fds_tmp, "Predios_SR")
        self.project_to_sr(FC_PREDIOS, pred_out, SR_OUT)
        self.compute_area_ha(pred_out)

        # 3) proyectar criterios presentes en cfg
        criterios_src = []
        for base in list_fcs_in_fds(FDS_CRITERIOS):
            if norm_name(base) not in cfg:
                continue
            src = os.path.join(FDS_CRITERIOS, base)
            dst = os.path.join(fds_tmp, "{}_SR".format(base))
            self.project_to_sr(src, dst, SR_OUT)
            criterios_src.append((base, dst))

        if not criterios_src:
            raise RuntimeError("No hay coincidencia entre tbl_config y las FC del FDS de criterios.")

        # 4) salida
        if not arcpy.Exists(FDS_SALIDA):
            arcpy.CreateFeatureDataset_management(out_gdb, os.path.basename(FDS_SALIDA), SR_OUT)
        out_fc = os.path.join(FDS_SALIDA, FC_SALIDA_NOM)
        if arcpy.Exists(out_fc):
            arcpy.Delete_management(out_fc)
        arcpy.CopyFeatures_management(pred_out, out_fc)

        # 5) procesar cada criterio
        for base, fc_path in criterios_src:
            k = norm_name(base)
            cfg_k = cfg.get(k, None)
            if not cfg_k:
                arcpy.AddWarning(u"Omitido: '{}' (sin fila en tbl_config)".format(base))
                continue

            metodo   = (cfg_k.get("metodo") or "minmax").lower()
            polarity = (cfg_k.get("polarity") or "benefit").lower()
            vmin     = cfg_k.get("vmin", 0.0)
            vmax     = cfg_k.get("vmax", 1.0)
            constv   = cfg_k.get("const")
            try:
                peso = float(cfg_k.get("peso", 0.0))
            except:
                peso = 0.0

            if not peso or peso == 0.0:
                arcpy.AddMessage(u"Omitido por peso=0: {}".format(base))
                continue

            cfg_i = {"metodo": metodo, "polarity": polarity, "vmin": vmin, "vmax": vmax, "const": constv}

            scores, covers = self.intersect_score(
                pred_out, fc_path, CAMPO_ID_PRED, CAMPO_VALOR, cfg_i
            )

            # campos de salida (Comp_ ahora es FIRMADO)
            f_score = arcpy.ValidateFieldName(u"Score_{}".format(base), out_fc)
            f_pct   = arcpy.ValidateFieldName(u"Pct_{}".format(base), out_fc)
            f_peso  = arcpy.ValidateFieldName(u"Peso_{}".format(base), out_fc)
            f_vmin  = arcpy.ValidateFieldName(u"Vmin_{}".format(base), out_fc)
            f_vmax  = arcpy.ValidateFieldName(u"Vmax_{}".format(base), out_fc)
            f_comp  = arcpy.ValidateFieldName(u"Comp_{}".format(base), out_fc)
            f_met   = arcpy.ValidateFieldName(u"Metodo_{}".format(base), out_fc)
            f_pol   = arcpy.ValidateFieldName(u"Polarity_{}".format(base), out_fc)

            for fn, tp in [(f_score,"DOUBLE"), (f_pct,"DOUBLE"), (f_peso,"DOUBLE"),
                           (f_vmin,"DOUBLE"), (f_vmax,"DOUBLE"), (f_comp,"DOUBLE"),
                           (f_met,"TEXT"), (f_pol,"TEXT")]:
                self.add_field_if_absent(out_fc, fn, tp)

            sign = -1.0 if polarity == "cost" else 1.0

            vmin_out = 0.0 if vmin is None else float(vmin)
            vmax_out = 1.0 if vmax is None else float(vmax)

            with arcpy.da.UpdateCursor(out_fc, [CAMPO_ID_PRED, f_score, f_pct, f_peso, f_vmin, f_vmax, f_comp, f_met, f_pol]) as cur:
                for pid, s, p, _, _, _, _, _, _ in cur:
                    sc   = scores.get(pid, 0.0)      # 0..1
                    pc   = covers.get(pid, 0.0)      # 0..1
                    comp = sc * peso * sign          # <<< AQUI: beneficio suma, costo resta
                    cur.updateRow([pid, sc, pc, peso, vmin_out, vmax_out, comp, metodo, polarity])

            arcpy.AddMessage(u"Incluido: {} (peso={}, metodo={}, pol={}, vmin={}, vmax={})"
                             .format(base, peso, metodo, polarity, vmin_out, vmax_out))

        # 6) totales y neto (nuevo)
        self.add_field_if_absent(out_fc, "Benef_Total", "DOUBLE")
        self.add_field_if_absent(out_fc, "Cost_Total",  "DOUBLE")
        self.add_field_if_absent(out_fc, "Indice_Neto", "DOUBLE")

        comp_fields = [f.name for f in arcpy.ListFields(out_fc) if f.name.upper().startswith("COMP_")]

        with arcpy.da.UpdateCursor(out_fc, comp_fields + ["Benef_Total", "Cost_Total", "Indice_Neto"]) as cur:
            n = len(comp_fields)
            for row in cur:
                comps = []
                for i in range(n):
                    try:
                        v = 0.0 if row[i] is None else float(row[i])
                    except:
                        v = 0.0
                    comps.append(v)
                benef = sum([c for c in comps if c > 0.0])
                cost  = sum([-c for c in comps if c < 0.0])  # valor absoluto
                neto  = benef - cost
                row[-3] = benef
                row[-2] = cost
                row[-1] = neto
                cur.updateRow(row)

        # 7) ranking por Indice_Neto
        self.add_field_if_absent(out_fc, "Rank_WLC", "LONG")
        rows = []
        with arcpy.da.SearchCursor(out_fc, ["OID@", "Indice_Neto"]) as cur:
            for oid, idx in cur:
                try:
                    v = 0.0 if idx in (None, "") else float(idx)
                except:
                    v = 0.0
                rows.append((oid, v))
        rows.sort(key=lambda r: r[1], reverse=True)
        rmap = {oid: i + 1 for i, (oid, _) in enumerate(rows)}
        with arcpy.da.UpdateCursor(out_fc, ["OID@", "Rank_WLC"]) as cur:
            for oid, _ in cur:
                cur.updateRow([oid, rmap.get(oid)])

        arcpy.AddMessage("Listo. Salida: {}".format(out_fc))

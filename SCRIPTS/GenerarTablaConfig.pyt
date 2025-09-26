# -*- coding: utf-8 -*-
import arcpy, os, sys, re, time, unicodedata
arcpy.env.overwriteOutput = True

# --------------- utils ---------------
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
    return str(b)

def strip_accents(s):
    s = to_text(s)
    return u"".join(c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn')

def norm_name(s):
    if not s:
        return u""
    s = strip_accents(to_text(s)).strip().upper()
    for suf in ["_AOI_WGS84","_WGS84","_AOI","_SR"]:
        if s.endswith(suf):
            s = s[:-len(suf)]
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s

def list_fcs_in_fds(fds_path):
    fcs = []
    if not fds_path or not arcpy.Exists(fds_path):
        return fcs
    for dp, dn, fns in arcpy.da.Walk(fds_path, datatype="FeatureClass"):
        for fn in fns:
            try:
                nm = arcpy.Describe(os.path.join(dp, fn)).baseName
            except:
                nm = os.path.basename(fn)
            # evitar capas de salida tipicas
            if nm.upper().startswith("PREDIOS_"):
                continue
            fcs.append(nm)
    return sorted(set(fcs))

def compute_minmax(fc_path, valor_field):
    vmin, vmax, n = None, None, 0
    if not arcpy.Exists(fc_path):
        return 0.0, 1.0, 0
    try:
        with arcpy.da.SearchCursor(fc_path, [valor_field]) as cur:
            for (v,) in cur:
                if v is None:
                    continue
                try:
                    fv = float(v)
                except:
                    continue
                n += 1
                vmin = fv if vmin is None else min(vmin, fv)
                vmax = fv if vmax is None else max(vmax, fv)
    except Exception as e:
        arcpy.AddWarning("No se pudo leer el campo {} en {}: {}".format(valor_field, fc_path, e))
        return 0.0, 1.0, 0
    if n == 0 or vmin is None or vmax is None or vmax == vmin:
        return 0.0, 1.0, n
    return float(vmin), float(vmax), int(n)

def ensure_table(gdb_path, table_name, overwrite_all):
    tbl_path = os.path.join(gdb_path, table_name)
    if arcpy.Exists(tbl_path) and overwrite_all:
        arcpy.Delete_management(tbl_path)
    if not arcpy.Exists(tbl_path):
        arcpy.CreateTable_management(gdb_path, table_name)
        # esquema
        arcpy.AddField_management(tbl_path, "Criterio", "TEXT", field_length=128)
        arcpy.AddField_management(tbl_path, "Criterio_norm", "TEXT", field_length=128)
        arcpy.AddField_management(tbl_path, "Metodo", "TEXT", field_length=16)     # minmax|max|binario|const
        arcpy.AddField_management(tbl_path, "Polarity", "TEXT", field_length=16)   # benefit|cost
        arcpy.AddField_management(tbl_path, "Vmin", "DOUBLE")
        arcpy.AddField_management(tbl_path, "Vmax", "DOUBLE")
        arcpy.AddField_management(tbl_path, "ConstVal", "DOUBLE")
        arcpy.AddField_management(tbl_path, "Peso", "DOUBLE")
        arcpy.AddField_management(tbl_path, "N", "LONG")
        arcpy.AddField_management(tbl_path, "FechaCalc", "TEXT", field_length=32)
        arcpy.AddField_management(tbl_path, "Fuente", "TEXT", field_length=32)
    return tbl_path

def read_existing_rows(tbl_path):
    prev = {}
    if not tbl_path or not arcpy.Exists(tbl_path):
        return prev
    fields = [f.name for f in arcpy.ListFields(tbl_path)]
    needed = ["Criterio","Criterio_norm","Metodo","Polarity","Vmin","Vmax","ConstVal","Peso"]
    for k in needed:
        if k not in fields:
            arcpy.AddWarning("Tabla existente sin campo {}. No se preserva.".format(k))
            return {}
    with arcpy.da.SearchCursor(tbl_path, needed) as cur:
        for C, Cn, M, P, Vmin, Vmax, CVal, W in cur:
            if not Cn:
                continue
            prev[to_text(Cn)] = {"Criterio": to_text(C),
                                 "Metodo": to_text(M) if M else "",
                                 "Polarity": to_text(P) if P else "",
                                 "Vmin": Vmin, "Vmax": Vmax,
                                 "ConstVal": CVal, "Peso": W}
    return prev

def delete_rows_for(tbl_path, crit_norm_set):
    if not crit_norm_set:
        return
    fld = "Criterio_norm"
    vals = list(crit_norm_set)
    step = 800
    i = 0
    while i < len(vals):
        chunk = vals[i:i+step]
        where = " OR ".join(["{} = '{}'".format(fld, v.replace("'", "''")) for v in chunk])
        with arcpy.da.UpdateCursor(tbl_path, [fld], where) as cur:
            for _ in cur:
                cur.deleteRow()
        i += step

def guess_defaults(fc_name_norm):
    n = fc_name_norm
    if re.search(r"(MINER|TITULO)", n):
        return ("binario", "cost")
    if re.search(r"(RAMSAR|SINAP|PNR|PROTEGID)", n):
        return ("binario", "benefit")
    if re.search(r"(LRE|ESTRUCTURA|CONECTIVI|BOSQUE|BST)", n):
        return ("minmax", "benefit")
    return ("minmax", "benefit")

# --------------- toolbox ---------------
class Toolbox(object):
    def __init__(self):
        self.label = "Generar tabla de configuracion"
        self.alias = "gen_tbl_config"
        self.tools = [GenerarTablaConfig]

class GenerarTablaConfig(object):
    def __init__(self):
        self.label = "Construir/actualizar tbl_config desde FDS"
        self.description = "Crea o actualiza una tabla con metodo, polarity, rangos y peso por criterio."
        self.canRunInBackground = True

    def getParameterInfo(self):
        p_gdb  = arcpy.Parameter(displayName="Geodatabase (GDB)", name="gdb",
                                 datatype="DEWorkspace", parameterType="Required", direction="Input")
        p_fds  = arcpy.Parameter(displayName="Feature Dataset de criterios", name="fds_criterios",
                                 datatype="DEFeatureDataset", parameterType="Required", direction="Input")
        p_fds.parameterDependencies = [p_gdb.name]
        p_val  = arcpy.Parameter(displayName="Campo 'Valor' en criterios", name="campo_valor",
                                 datatype="GPString", parameterType="Required", direction="Input")
        p_val.value = "Valor"
        p_tbl  = arcpy.Parameter(displayName="Nombre de tabla (en GDB)", name="nombre_tabla",
                                 datatype="GPString", parameterType="Required", direction="Input")
        p_tbl.value = "tbl_config"
        p_over = arcpy.Parameter(displayName="Sobrescribir tabla completa", name="overwrite_all",
                                 datatype="GPBoolean", parameterType="Optional", direction="Input")
        p_over.value = False
        p_pres = arcpy.Parameter(displayName="Preservar edicion previa", name="preservar_manual",
                                 datatype="GPBoolean", parameterType="Optional", direction="Input")
        p_pres.value = True
        p_norm = arcpy.Parameter(displayName="Normalizar pesos (suma 1)", name="normalizar_pesos",
                                 datatype="GPBoolean", parameterType="Optional", direction="Input")
        p_norm.value = True
        p_defw = arcpy.Parameter(displayName="Peso por defecto (0-1)", name="peso_default",
                                 datatype="GPDouble", parameterType="Optional", direction="Input")
        p_defw.value = 0.0
        p_inf = arcpy.Parameter(displayName="Inferir metodo/polarity por nombre", name="inferir_metodos",
                                datatype="GPBoolean", parameterType="Optional", direction="Input")
        p_inf.value = True
        p_list = arcpy.Parameter(displayName="Lista manual de FC (coma)", name="fc_lista",
                                 datatype="GPString", parameterType="Optional", direction="Input")
        return [p_gdb, p_fds, p_val, p_tbl, p_over, p_pres, p_norm, p_defw, p_inf, p_list]

    def updateMessages(self, params):
        if params[0].valueAsText and not arcpy.Exists(params[0].valueAsText):
            params[0].setErrorMessage("La GDB no existe.")
        if params[1].valueAsText and not arcpy.Exists(params[1].valueAsText):
            params[1].setErrorMessage("El FDS no existe.")

    def execute(self, params, messages):
        GDB           = params[0].valueAsText
        FDS           = params[1].valueAsText
        CAMPO_VALOR   = params[2].valueAsText
        NOMBRE_TABLA  = params[3].valueAsText
        OVERWRITE_ALL = bool(params[4].value)
        PRESERVAR     = bool(params[5].value)
        NORMALIZAR    = bool(params[6].value)
        PESO_DEFAULT  = float(params[7].value) if params[7].value is not None else 0.0
        INFERIR       = bool(params[8].value)
        FC_LISTA      = params[9].valueAsText

        if FC_LISTA:
            fcs = [x.strip() for x in FC_LISTA.split(",") if x.strip()]
        else:
            fcs = list_fcs_in_fds(FDS)
        if not fcs:
            raise RuntimeError("No se encontraron FC en el FDS.")

        tbl_path = ensure_table(GDB, NOMBRE_TABLA, OVERWRITE_ALL)
        prev = read_existing_rows(tbl_path) if (PRESERVAR and not OVERWRITE_ALL) else {}
        fecha = time.strftime("%Y-%m-%d %H:%M:%S")

        rows = []
        pesos_tmp = []
        for fc in fcs:
            fc_path = os.path.join(FDS, fc)
            if not arcpy.Exists(fc_path):
                arcpy.AddWarning("FC no encontrada: {} (omitida)".format(fc))
                continue

            k = norm_name(fc)
            vmin, vmax, n = compute_minmax(fc_path, CAMPO_VALOR)
            metodo, pol = ("minmax", "benefit")
            if INFERIR:
                metodo, pol = guess_defaults(k)
            peso = PESO_DEFAULT
            constv = None

            if k in prev:
                if prev[k].get("Metodo"):   metodo = to_text(prev[k]["Metodo"])
                if prev[k].get("Polarity"): pol    = to_text(prev[k]["Polarity"])
                if prev[k].get("Vmin") is not None:    vmin   = float(prev[k]["Vmin"])
                if prev[k].get("Vmax") is not None:    vmax   = float(prev[k]["Vmax"])
                if prev[k].get("ConstVal") is not None: constv = float(prev[k]["ConstVal"])
                if prev[k].get("Peso") is not None:    peso   = float(prev[k]["Peso"])

            if constv is None and metodo.lower() == "const":
                constv = 1.0

            rows.append([fc, k, metodo, pol, vmin, vmax, constv, peso, n, fecha, ("previo" if k in prev else "auto")])
            pesos_tmp.append(peso)

        if NORMALIZAR:
            s = sum(pesos_tmp)
            if s > 0:
                for i in range(len(rows)):
                    rows[i][7] = rows[i][7] / s

        if not OVERWRITE_ALL:
            delete_rows_for(tbl_path, set(norm_name(fc) for fc in fcs))

        fields = ["Criterio","Criterio_norm","Metodo","Polarity","Vmin","Vmax","ConstVal","Peso","N","FechaCalc","Fuente"]
        with arcpy.da.InsertCursor(tbl_path, fields) as icur:
            for r in rows:
                icur.insertRow(r)
                arcpy.AddMessage("OK {} -> metodo={} polarity={} vmin={} vmax={} const={} peso={} n={}".format(
                    r[0], r[2], r[3], r[4], r[5], r[6], r[7], r[8]
                ))
        arcpy.AddMessage("Tabla lista: {}".format(tbl_path))

# -*- coding: utf-8 -*-
# Compatibilidad: ArcGIS Desktop 10.8 (Python 2.7)
import arcpy, os, sys, unicodedata, re, time

arcpy.env.overwriteOutput = True

# ---------------- utilidades ----------------
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

def gdb_name_from_path(gdb_path):
    try:
        return os.path.basename(gdb_path.rstrip("\\/"))
    except:
        return to_text(gdb_path)

# ---------------- listar FC en FDS ----------------
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
            # evitar capas tipicas de salida
            if nm and nm.upper().startswith("PREDIOS_"):
                continue
            fcs.append(nm)
    return sorted(set([x for x in fcs if x]))

# ---------------- verificacion campo Valor ----------------
_NUMERIC_TYPES = set(["Double","Single","Integer","SmallInteger"])

def get_valor_field_if_numeric(fc_path, field_name="Valor"):
    """Retorna el nombre real del campo si existe y es numerico; si no, None."""
    try:
        for f in arcpy.ListFields(fc_path):
            if f.name.lower() == field_name.lower():
                return f.name if f.type in _NUMERIC_TYPES else None
    except Exception as e:
        arcpy.AddWarning(u"No se pudo listar campos en '{}': {}".format(fc_path, to_text(e)))
        return None
    return None

# ---------------- estadisticos e inferencia ----------------
def stats_and_method_by_threshold(fc_path, valor_field):
    """
    Calcula vmin, vmax, n y define metodo por umbral:
      - si vmax <= 1.0 => Metodo='binario' (Vmin=0, Vmax=1)
      - si vmax  > 1.0 => Metodo='minmax' (usa Vmin/Vmax; si vmin==vmax => 0..1)
    Retorna (metodo, vmin, vmax, n)
    """
    n = 0
    vmin, vmax = None, None
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
        arcpy.AddWarning(u"No se pudo leer '{}' en '{}': {}".format(valor_field, fc_path, to_text(e)))
        return ("minmax", 0.0, 1.0, 0)

    if n == 0:
        return ("minmax", 0.0, 1.0, 0)

    if vmax is not None and vmax <= 1.0:
        return ("binario", 0.0, 1.0, n)

    if vmin == vmax:
        return ("minmax", 0.0, 1.0, n)

    return ("minmax", float(vmin), float(vmax), n)

# ---------------- crear/recrear tabla ----------------
def recreate_config_table(output_gdb, table_name):
    # validar que el destino no sea un feature dataset
    try:
        desc = arcpy.Describe(output_gdb)
        if getattr(desc, "dataType", "").lower() == "featuredataset":
            raise arcpy.ExecuteError(u"El destino seleccionado es un Feature Dataset. "
                                     u"Las tablas deben crearse en la GDB raiz.")
    except:
        pass

    tbl_path = os.path.join(output_gdb, table_name)
    if arcpy.Exists(tbl_path):
        arcpy.Delete_management(tbl_path)
    arcpy.CreateTable_management(output_gdb, table_name)

    # esquema (sin ConstVal) + alias/descripcion sin tildes
    arcpy.AddField_management(tbl_path, "Criterio",      "TEXT",   field_length=128, field_alias="Criterio (nombre legible)")
    arcpy.AddField_management(tbl_path, "Criterio_norm", "TEXT",   field_length=128, field_alias="Criterio_norm (clave normalizada)")
    arcpy.AddField_management(tbl_path, "Metodo",        "TEXT",   field_length=16,  field_alias="Metodo (binario|minmax|max)")
    arcpy.AddField_management(tbl_path, "Polarity",      "TEXT",   field_length=16,  field_alias="Polarity (benefit|cost)")
    arcpy.AddField_management(tbl_path, "Vmin",          "DOUBLE",                    field_alias="Vmin (minimo observado o 0)")
    arcpy.AddField_management(tbl_path, "Vmax",          "DOUBLE",                    field_alias="Vmax (maximo observado o 1)")
    arcpy.AddField_management(tbl_path, "Peso",          "DOUBLE",                    field_alias="Peso (ponderacion del criterio)")
    arcpy.AddField_management(tbl_path, "N",             "LONG",                      field_alias="N (conteo de valores validos)")
    arcpy.AddField_management(tbl_path, "FechaCalc",     "TEXT",   field_length=32,   field_alias="FechaCalc (timestamp de calculo)")
    arcpy.AddField_management(tbl_path, "Fuente",        "TEXT",   field_length=64,   field_alias="Fuente (nombre de la GDB de entrada)")
    return tbl_path

# ---------------- toolbox ----------------
class Toolbox(object):
    def __init__(self):
        self.label = "Generar tabla de configuracion (auto desde FDS)"
        self.alias = "gen_tbl_config_auto"
        self.tools = [GenerarTablaConfigAuto]

class GenerarTablaConfigAuto(object):
    def __init__(self):
        self.label = "Crear/llenar tbl_config leyendo 'Valor' de cada FC"
        self.description = (
            "ArcGIS Desktop 10.8 (Py2.7). Recorre todo el FDS y crea la tabla 'tbl_config' en la GDB de salida. "
            "Nota importante: cada Feature Class de entrada debe tener un campo numerico llamado 'Valor' (obligatorio). "
            "Inferencia de Metodo por umbral del maximo: vmax <= 1 => binario (Vmin=0,Vmax=1); vmax > 1 => minmax "
            "(usa Vmin/Vmax observados; si vmin==vmax => 0..1). Las FC sin 'Valor' numerico se omiten con AddWarning. "
            "El campo 'Fuente' guarda el nombre de la GDB de entrada."
        )
        self.canRunInBackground = True

    def getParameterInfo(self):
        p_in_gdb = arcpy.Parameter(
            displayName="GDB de entrada (contiene el FDS)",
            name="in_gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input"
        )

        p_fds = arcpy.Parameter(
            displayName="Feature Dataset de criterios (entrada) - Nota: es obligatorio que cada FC tenga un campo 'Valor' numerico",
            name="fds_criterios",
            datatype="DEFeatureDataset",
            parameterType="Required",
            direction="Input"
        )
        p_fds.parameterDependencies = [p_in_gdb.name]

        p_out_gdb = arcpy.Parameter(
            displayName="GDB de salida (crear tabla aqui)",
            name="out_gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input"
        )

        p_tbl = arcpy.Parameter(
            displayName="Nombre de la tabla (en la GDB de salida)",
            name="nombre_tabla",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        p_tbl.value = "tbl_config"

        return [p_in_gdb, p_fds, p_out_gdb, p_tbl]

    def updateMessages(self, params):
        if params[0].valueAsText and not arcpy.Exists(params[0].valueAsText):
            params[0].setErrorMessage("La GDB de entrada no existe.")
        if params[1].valueAsText and not arcpy.Exists(params[1].valueAsText):
            params[1].setErrorMessage("El FDS de entrada no existe.")
        if params[2].valueAsText and not arcpy.Exists(params[2].valueAsText):
            params[2].setErrorMessage("La GDB de salida no existe.")
        if params[2].valueAsText and arcpy.Exists(params[2].valueAsText):
            try:
                dt = getattr(arcpy.Describe(params[2].valueAsText), "dataType", "")
                if str(dt).lower() == "featuredataset":
                    params[2].setErrorMessage("La salida no puede ser un Feature Dataset. Seleccione la GDB raiz.")
            except:
                pass

    def execute(self, params, messages):
        IN_GDB     = params[0].valueAsText
        FDS        = params[1].valueAsText
        OUT_GDB    = params[2].valueAsText
        TABLE_NAME = params[3].valueAsText

        fuente_gdb = gdb_name_from_path(IN_GDB)

        # 1) recrear tabla en la GDB de salida
        tbl_path = recreate_config_table(OUT_GDB, TABLE_NAME)

        # 2) recorrer todo el FDS
        fcs = list_fcs_in_fds(FDS)
        if not fcs:
            raise RuntimeError("No se encontraron Feature Class en el FDS.")

        fecha = time.strftime("%Y-%m-%d %H:%M:%S")
        fields = ["Criterio","Criterio_norm","Metodo","Polarity","Vmin","Vmax","Peso","N","FechaCalc","Fuente"]

        inserted = 0
        omitted = []

        with arcpy.da.InsertCursor(tbl_path, fields) as icur:
            for fc in fcs:
                fc_path = os.path.join(FDS, fc)
                if not arcpy.Exists(fc_path):
                    arcpy.AddWarning(u"FC no encontrada: {} (omitida)".format(fc))
                    omitted.append(fc)
                    continue

                valor_field = get_valor_field_if_numeric(fc_path, "Valor")
                if not valor_field:
                    arcpy.AddWarning(u"La FC '{}' no tiene campo 'Valor' numerico. Se omite.".format(fc))
                    omitted.append(fc)
                    continue

                metodo, vmin, vmax, n = stats_and_method_by_threshold(fc_path, valor_field)

                pol  = "benefit"
                peso = 1.0
                k    = norm_name(fc)

                icur.insertRow([fc, k, metodo, pol, vmin, vmax, peso, n, fecha, fuente_gdb])
                inserted += 1
                arcpy.AddMessage(u"OK {:<30} -> metodo={} vmin={} vmax={} n={}".format(fc, metodo, vmin, vmax, n))

        arcpy.AddMessage(u"Tabla creada en: {}".format(tbl_path))
        arcpy.AddMessage(u"Filas insertadas: {}".format(inserted))
        if omitted:
            arcpy.AddWarning(u"FC omitidas (sin 'Valor' numerico o no encontradas): {}".format(", ".join(omitted)))
        arcpy.AddMessage(u"Nota: revise y edite manualmente 'Polarity' y 'Peso' segun su modelo.")

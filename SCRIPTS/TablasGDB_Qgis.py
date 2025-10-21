# QGIS 3.x - PyQGIS
# -*- coding: utf-8 -*-
# Ejecuta en la Consola Python de QGIS o como script en el Editor de QGIS

from qgis.core import (
    QgsVectorLayer, QgsFields, QgsField, QgsFeature, QgsWkbTypes,
    QgsVectorFileWriter, QgsCoordinateReferenceSystem
)
from PyQt5.QtCore import QVariant
import os, datetime

# =============== CONFIG ===============
EXCEL_PATH     = r"D:\XLS_Temporal\MuestreoFloraResultados_ODL_BHM.xlsx"
SHEET_DATOS    = "Datos"        
SHEET_DICT     = "Diccionario"  

# Salida como DBF (sin geometría), recomendado
OUT_FOLDER     = r"D:\XLS_Temporal\salidas_dbf"
OUT_DBF_NAME   = "Seg_EspSembradaTB_CLEAN.dbf"

# (Opcional) Salida directa a FileGDB si tienes el driver FileGDB habilitado:
USE_FILEGDB    = False
OUT_GDB        = r"D:\XLS_Temporal\GDB_ANLA_sandbox.gdb"
OUT_TABLE_NAME = "Seg_EspSembradaTB_CLEAN323232"
# =====================================

# Mapa de tipos ANLA -> PyQGIS/OGR
TYPE_MAP = {
    "single":      ("REAL",   7, 6),   # (tipo, precision, escala para redondeo)
    "double":      ("REAL",  15, 8),
    "smallinteger":("INT",   None, None),
    "longinteger": ("INT",   None, None),
    "string":      ("TEXT",  None, None),  # longitud viene en TAMAÑO
    "date":        ("DATE",  None, None)
}

# Defaults si no hay diccionario (opcional)
TYPE_MAP_DEFAULTS = {
    # "DEN_SIEMB": ("REAL", 7, 6),
    # "MANT_PER":  ("INT", None, None),
}

def log(msg):
    print(msg)

def ensure_folder(path):
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

def excel_layer(path, sheet):
    """Crea un QgsVectorLayer apuntando a una hoja Excel (GDAL/OGR)."""
    uri = r'{}|layername={}'.format(path, sheet)  # ruta.xlsx|layername=Hoja
    vl = QgsVectorLayer(uri, sheet, "ogr")
    if not vl.isValid():
        raise RuntimeError("No se pudo abrir la hoja de Excel: {}".format(uri))
    return vl

def read_schema_from_dict(dict_vl):
    """Lee CAMPO, TIPO_DE_DATO, TAMAÑO desde la hoja Diccionario ya cargada como capa."""
    idx_campo = dict_vl.fields().indexFromName("CAMPO")
    idx_tipo  = dict_vl.fields().indexFromName("TIPO_DE_DATO")
    idx_tam   = dict_vl.fields().indexFromName("TAMAÑO")
    if min(idx_campo, idx_tipo, idx_tam) < 0:
        raise RuntimeError("La hoja 'Diccionario' debe contener las columnas CAMPO, TIPO_DE_DATO y TAMAÑO.")

    schema = {}
    for f in dict_vl.getFeatures():
        campo = (f[idx_campo] or "").strip()
        tipo  = (f[idx_tipo] or "").strip().lower()
        tam   = f[idx_tam]
        if not campo:
            continue
        if tipo not in TYPE_MAP:
            raise ValueError("Tipo no soportado en diccionario: {}".format(tipo))
        _kind, prec_def, scale_def = TYPE_MAP[tipo]
        length = None
        if _kind == "TEXT":
            try:
                length = int(tam) if tam not in (None, "", " ") else 255
            except:
                length = 255
        schema[campo] = {
            "kind": _kind,     # TEXT | INT | REAL | DATE
            "precision": prec_def,
            "scale": scale_def,
            "length": length
        }
    return schema

def merge_defaults(schema, defaults):
    sch = dict(schema) if schema else {}
    for k, (tp, pr, sc) in defaults.items():
        if k not in sch:
            sch[k] = {
                "kind": ("TEXT" if tp=="TEXT" else ("REAL" if tp in ("DOUBLE","FLOAT","REAL") else "INT")),
                "precision": pr, "scale": sc, "length": (255 if tp=="TEXT" else None)
            }
    return sch

# --- helpers de coerción ---
def to_float(v):
    if v in (None, ""): return None
    try:
        if isinstance(v, str):
            v = v.replace(",", ".")
        return float(v)
    except:
        return None

def to_int(v):
    if v in (None, ""): return None
    try:
        if isinstance(v, str):
            v = v.replace(",", ".")
        return int(round(float(v)))
    except:
        return None

def to_text(v, n):
    if v in (None, ""): return None
    s = str(v)
    return s[:n] if n else s

def to_date(v):
    # QGIS/OGR Date guarda solo fecha (sin hora).
    if v in (None, ""): return None
    if isinstance(v, datetime.date) and not isinstance(v, datetime.datetime):
        return v
    if isinstance(v, datetime.datetime):
        return v.date()
    fmts = ["%Y-%m-%d","%d/%m/%Y","%d-%m-%Y","%Y/%m/%d","%d/%m/%Y %H:%M:%S"]
    for f in fmts:
        try:
            return datetime.datetime.strptime(str(v), f).date()
        except:
            pass
    return None

def round_to(v, n):
    if v is None: return None
    try:
        return round(float(v), n)
    except:
        return None

def build_field_definitions(schema):
    fields = QgsFields()
    for name, spec in schema.items():
        kind = spec["kind"]
        if kind == "TEXT":
            fields.append(QgsField(name, QVariant.String, spec.get("length", 255)))
        elif kind == "INT":
            fields.append(QgsField(name, QVariant.Int))
        elif kind == "REAL":
            fields.append(QgsField(name, QVariant.Double))
        elif kind == "DATE":
            fields.append(QgsField(name, QVariant.Date))
        else:
            raise ValueError("Tipo no soportado: {}".format(kind))
    return fields

def build_row_converter(in_fields, schema):
    out_fields = list(schema.keys())
    index = {f.name(): i for i, f in enumerate(in_fields)}

    funcs = []
    for f in out_fields:
        spec = schema[f]
        kind = spec["kind"]
        if f not in index:
            funcs.append(lambda feat, _f=f: None)
            continue
        i = index[f]
        if kind == "TEXT":
            L = spec.get("length", 255)
            funcs.append(lambda feat, ii=i, LL=L: to_text(feat[ii], LL))
        elif kind == "REAL":
            decs = spec.get("scale", 6) or 6
            funcs.append(lambda feat, ii=i, dd=decs: round_to(to_float(feat[ii]), dd))
        elif kind == "INT":
            funcs.append(lambda feat, ii=i: to_int(feat[ii]))
        elif kind == "DATE":
            funcs.append(lambda feat, ii=i: to_date(feat[ii]))
        else:
            funcs.append(lambda feat, ii=i: feat[ii])

    def convert(feat):
        return [fn(feat) for fn in funcs]

    return convert, out_fields

def main_qgis():
    ensure_folder(OUT_FOLDER)

    log("=== Limpiador XLSX → Tabla (.dbf o FileGDB) con tipos ANLA ===")
    log("Excel: {}".format(EXCEL_PATH))

    # 1) Cargar hojas de Excel como capas OGR
    vl_datos = excel_layer(EXCEL_PATH, SHEET_DATOS)

    # Diccionario (opcional)
    schema = {}
    try:
        vl_dict  = excel_layer(EXCEL_PATH, SHEET_DICT)
        schema   = read_schema_from_dict(vl_dict)
        log("Diccionario leído: {} campos".format(len(schema)))
    except Exception as e:
        log("(Aviso) No se pudo leer 'Diccionario': {}. Usaré TYPE_MAP_DEFAULTS y/o inferencia TEXT(255).".format(e))

    # 2) Defaults y/o inferencia
    schema = merge_defaults(schema, TYPE_MAP_DEFAULTS)
    if not schema:
        schema = {f.name(): {"kind":"TEXT", "length":255, "precision":None, "scale":None}
                  for f in vl_datos.fields() if f.name().lower() not in ("fid","ogc_fid")}
        log("(Aviso) Sin diccionario: todos los campos serán TEXT(255).")

    # 3) Crear estructura de salida (sin geometría)
    out_fields = build_field_definitions(schema)

    # 4) Preparar escritor (DBF o FileGDB)
    dest_path = None
    if USE_FILEGDB:
        # Requiere driver FileGDB en tu GDAL/QGIS
        from qgis.core import QgsVectorFileWriter
        options = QgsVectorFileWriter.SaveVectorOptions()
        options.driverName = "FileGDB"
        options.layerName  = OUT_TABLE_NAME
        options.fileEncoding = "UTF-8"
        # En FileGDB, el fileName es la .gdb (se crea capa interna)
        writer = QgsVectorFileWriter.create(
            OUT_GDB, out_fields, QgsWkbTypes.NoGeometry,
            QgsCoordinateReferenceSystem("EPSG:4326"),
            options
        )
        if writer.hasError() != QgsVectorFileWriter.NoError:
            raise RuntimeError("Error creando tabla en FileGDB: {}".format(writer.errorMessage()))
        dest_path = os.path.join(OUT_GDB, OUT_TABLE_NAME)
    else:
        # ESRI Shapefile (tabla .dbf sin geometría)
        dbf_path = os.path.join(OUT_FOLDER, OUT_DBF_NAME)
        options = QgsVectorFileWriter.SaveVectorOptions()
        options.driverName   = "ESRI Shapefile"
        options.fileEncoding = "UTF-8"
        writer = QgsVectorFileWriter.create(
            dbf_path, out_fields, QgsWkbTypes.NoGeometry,
            QgsCoordinateReferenceSystem("EPSG:4326"),
            options
        )
        if writer.hasError() != QgsVectorFileWriter.NoError:
            raise RuntimeError("Error creando DBF: {}".format(writer.errorMessage()))
        dest_path = dbf_path

    # 5) Insertar filas con conversión/coerción
    convert, out_names = build_row_converter(vl_datos.fields(), schema)
    total = 0
    for feat in vl_datos.getFeatures():
        vals = convert(feat)
        newf = QgsFeature()
        newf.setFields(out_fields)
        newf.initAttributes(out_fields.count())
        for idx, _name in enumerate(out_names):
            newf.setAttribute(idx, vals[idx])
        writer.addFeature(newf)
        total += 1

    del writer  # cierra/volca

    log("Salida: {} ({} filas insertadas)".format(dest_path, total))
    return dest_path

# Ejecuta:
if __name__ == "__main__":
    path = main_qgis()
    log("OK: {}".format(path))

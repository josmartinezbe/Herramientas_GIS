# -*- coding: utf-8 -*-
import arcpy, os, datetime

# ========= PARÁMETROS DEL USUARIO =========
# Ejecutar dentro de ArcMap (usa el MXD abierto)
mxd_path       = "D:\HERRAMIENTAS_GIS\MXD\250828_Frente_Mantenimiento_4_ODL2025.mxd"  # si prefieres un .mxd específico, pon su ruta

shp_path       = r"D:\HERRAMIENTAS_GIS\SHP\250828_ManteODL.shp"
lyr_template   = r"D:\HERRAMIENTAS_GIS\LYR\250828_ManteODL.lyr"
output_folder  = r"D:\HERRAMIENTAS_GIS\MAPAS"

layer_name_in_toc = "250828_ManteODL"  # cómo quieres que se llame la capa en el TOC
campo_rangos      = "ID_PG"       # << se usa para los filtros/rangos
campo_label       = "ID_Rest_25"       # << el .lyr debería etiquetar con este campo

# (Opcional) filtrar por EXPEDIENTE. Déjalo [] para no filtrar.
expedientes = []  # p.ej.: ["LAM2965"]

# Rangos: (min_exclusivo, max_inclusivo, alias_salida)
rangos = [
    (16, 45,  "16_45"),
    (46, 75,  "46_75"),
    (76, 105, "76_105"),
    (106,135, "106_135")
]

dpi_export = 300
# ========= FIN PARÁMETROS =========

# Utilidades
def ensure_folder(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def get_or_add_layer(df, shp, desired_name):
    # Buscar capa existente por dataSource o por nombre
    for lyr in arcpy.mapping.ListLayers(mxd, "", df):
        try:
            if lyr.isFeatureLayer and (lyr.dataSource.lower() == shp.lower() or lyr.name == desired_name):
                return lyr
        except:
            pass
    # Agregar si no existe
    new_lyr = arcpy.mapping.Layer(shp)
    arcpy.mapping.AddLayer(df, new_lyr, "TOP")
    # Renombrar en TOC
    new_lyr.name = desired_name
    return new_lyr

# Preparación
ensure_folder(output_folder)
mxd = arcpy.mapping.MapDocument(mxd_path)
df  = arcpy.mapping.ListDataFrames(mxd)[0]   # usa el primer dataframe

# Asegurar capa en el MXD
lyr = get_or_add_layer(df, shp_path, layer_name_in_toc)

# Cargar plantilla (labels: placement regular, offset horizontal 30 pt, bullet leader)
tmpl = arcpy.mapping.Layer(lyr_template)
arcpy.mapping.UpdateLayer(df, lyr, tmpl, True)

# Encender etiquetas por si la plantilla no lo dejó activo
try:
    lyr.showLabels = True
except:
    pass

# Construir cláusula base (por EXPEDIENTE si aplica)
if expedientes:
    in_vals = ", ".join("'{}'".format(x) for x in expedientes)
    where_base = "EXPEDIENTE IN ({})".format(in_vals)
else:
    where_base = "1=1"

# Exportar un PDF por rango
for minimo, maximo, etiqueta in rangos:
    where_rango = "{wb} AND {campo} > {minv} AND {campo} <= {maxv}".format(
        wb=where_base, campo=campo_rangos, minv=minimo, maxv=maximo
    )
    lyr.definitionQuery = where_rango

    # Reaplicar plantilla (asegura el estilo de etiquetas en cada iteración)
    arcpy.mapping.UpdateLayer(df, lyr, tmpl, True)
    try:
        lyr.showLabels = True
    except:
        pass

    # Contar entidades del rango (para el nombre del archivo)
    tmp_view = "view_{}".format(etiqueta)
    arcpy.MakeFeatureLayer_management(lyr.dataSource, tmp_view, where_rango)
    recs = int(arcpy.GetCount_management(tmp_view).getOutput(0))
    arcpy.Delete_management(tmp_view)

    stamp  = datetime.datetime.now().strftime("%Y%m%d")
    outpdf = os.path.join(output_folder, u"Mapa_{0}_{1}regs_{2}.pdf".format(etiqueta, recs, stamp))
    arcpy.mapping.ExportToPDF(mxd, outpdf, df, resolution=dpi_export)
    print(u"Exportado: {} ({} registros)".format(outpdf, recs))

del mxd
print("Hecho.")

# -*- coding: utf-8 -*-
import arcpy, os, datetime
####################################################################################################
# ===== PARÁMETROS (reemplaza mxd_path con tu ruta) =====
mxd_path       = r"D:\HERRAMIENTAS_GIS\MXD\250828_Frente_Mantenimiento_4_ODL2025.mxd"      # <-- tu MXD
shp_path       = r"D:\HERRAMIENTAS_GIS\SHP\250828_ManteODL.shp"
lyr_template   = r"D:\HERRAMIENTAS_GIS\LYR\250828_ManteODL.lyr"
output_folder  = r"D:\HERRAMIENTAS_GIS\MAPAS"

layer_name_in_toc = "250828_ManteODL"
campo_rangos      = "ID_PG"         # filtro por rangos
# etiquetas: las toma del .lyr (ID_Rest_25)
expedientes       = ["LAM2965"]     # [] para no filtrar por expediente

rangos = [(16,45,"16_45"), (46,75,"46_75"), (76,105,"76_105"), (106,135,"106_135")]
dpi_export = 300
# ========================================================

def ensure_folder(p):
    if not os.path.isdir(p): os.makedirs(p)

def base_where(expeds):
    if not expeds: return "1=1"
    vals = ", ".join("'{}'".format(x.replace("'", "''")) for x in expeds)
    return "EXPEDIENTE IN ({})".format(vals)

def get_or_add_layer(df, shp, desired_name):
    for lyr in arcpy.mapping.ListLayers(mxd, "", df):
        try:
            if lyr.isFeatureLayer and (lyr.dataSource.lower()==shp.lower() or lyr.name==desired_name):
                return lyr
        except: pass
    new_lyr = arcpy.mapping.Layer(shp)
    arcpy.mapping.AddLayer(df, new_lyr, "TOP")
    new_lyr = arcpy.mapping.ListLayers(mxd, new_lyr.name, df)[0]
    new_lyr.name = desired_name
    return new_lyr

ensure_folder(output_folder)

mxd = arcpy.mapping.MapDocument(mxd_path)
df  = arcpy.mapping.ListDataFrames(mxd)[0]  # primer dataframe, sin depender de "Layers"

# agregar capa si hace falta
lyr = get_or_add_layer(df, shp_path, layer_name_in_toc)

# aplicar plantilla (.lyr) con labels externos (offset 30 pt + bullet leader)
tmpl = arcpy.mapping.Layer(lyr_template)
arcpy.mapping.UpdateLayer(df, lyr, tmpl, True)
try: lyr.showLabels = True
except: pass

# zoom inicial al dataset completo (deja margen para leaders)
if lyr.isFeatureLayer and lyr.getExtent():
    df.extent = lyr.getExtent()
    df.scale  = df.scale * 1.1

where_base = base_where(expedientes)

for minimo, maximo, alias in rangos:
    where_rango = "{wb} AND {campo} > {minv} AND {campo} <= {maxv}".format(
        wb=where_base, campo=campo_rangos, minv=minimo, maxv=maximo
    )
    lyr.definitionQuery = where_rango
    arcpy.mapping.UpdateLayer(df, lyr, tmpl, True)
    try: lyr.showLabels = True
    except: pass

    # zoom al subconjunto si hay entidades (para que se vean bien los leaders)
    tmp_view = "view_{}".format(alias)
    arcpy.MakeFeatureLayer_management(lyr.dataSource, tmp_view, where_rango)
    recs = int(arcpy.GetCount_management(tmp_view).getOutput(0))
    if recs > 0:
        desc = arcpy.Describe(tmp_view)
        if hasattr(desc, "extent") and desc.extent:
            df.extent = desc.extent
            df.scale  = df.scale * 1.2
    arcpy.Delete_management(tmp_view)

    out_pdf = os.path.join(
        output_folder,
        u"Mapa_{0}_{1}regs_{2}.pdf".format(alias, recs, datetime.datetime.now().strftime("%Y%m%d"))
    )
    arcpy.mapping.ExportToPDF(mxd, out_pdf, df, resolution=dpi_export)
    print(u"Exportado: {} ({} registros)".format(out_pdf, recs))

mxd.save()   # opcional
del mxd
print("Listo.")


############################################################ -*- coding: utf-8 -*-
import arcpy, os, datetime, math

# ========= PARÁMETROS =========
mxd_path       = r"D:\HERRAMIENTAS_GIS\MXD\250828_Frente_Mantenimiento_4_ODL2025.mxd"      # <--- pon aquí tu MXD
shp_path       = r"D:\HERRAMIENTAS_GIS\SHP\250828_ManteODL.shp"
lyr_template   = r"D:\HERRAMIENTAS_GIS\LYR\250828_ManteODL.lyr"
output_folder  = r"D:\HERRAMIENTAS_GIS\MAPAS"

layer_name_in_toc = "250828_ManteODL"
campo_rangos      = "ID_PG"        # filtro por rangos (tu requisito)
# etiquetas: vienen del .lyr (debe rotular ID_Rest_25 con offset horizontal + bullet leader)

# Filtro base por expediente (vacío si no aplica)
expedientes = ["LAM2965"]

# Rangos (abierto por abajo, cerrado por arriba)
rangos = [
    (16, 45,  "16_45"),
    (46, 75,  "46_75"),
    (76, 105, "76_105"),
    (106,135, "106_135")
]

dpi_export   = 300
margen_ext   = 1.15   # expande un poco la extensión (para leaders)
# Lista de escalas "cerradas" permitidas (ajústala a tu cartografía)
escalas_cerradas = [
    2000, 2500, 3000, 4000, 5000, 6000, 7500, 10000, 12500, 15000, 20000,
    25000, 30000, 40000, 50000, 60000, 75000, 100000, 125000, 150000, 200000
]
# =============================

def ensure_folder(p):
    if not os.path.isdir(p):
        os.makedirs(p)

def base_where(expeds):
    if not expeds: return "1=1"
    vals = ", ".join("'{}'".format(x.replace("'", "''")) for x in expeds)
    return "EXPEDIENTE IN ({})".format(vals)

def get_or_add_layer(df, shp, desired_name):
    for lyr in arcpy.mapping.ListLayers(mxd, "", df):
        try:
            if lyr.isFeatureLayer and (lyr.dataSource.lower()==shp.lower() or lyr.name==desired_name):
                return lyr
        except:
            pass
    new_lyr = arcpy.mapping.Layer(shp)
    arcpy.mapping.AddLayer(df, new_lyr, "TOP")
    new_lyr = arcpy.mapping.ListLayers(mxd, new_lyr.name, df)[0]
    new_lyr.name = desired_name
    return new_lyr

def extent_expand(ext, factor):
    """Expande una extensión (arcpy.Extent) por un factor."""
    cx = (ext.XMin + ext.XMax) / 2.0
    cy = (ext.YMin + ext.YMax) / 2.0
    w  = (ext.XMax - ext.XMin) * factor
    h  = (ext.YMax - ext.YMin) * factor
    return arcpy.Extent(cx - w/2.0, cy - h/2.0, cx + w/2.0, cy + h/2.0)

def snap_escala(valor, escalas):
    """Devuelve la escala 'cerrada' más pequeña que sea >= valor.
       Si todas son menores, devuelve la mayor de la lista."""
    for s in escalas:
        if s >= valor:
            return s
    return escalas[-1]

# --- Preparación ---
ensure_folder(output_folder)

mxd = arcpy.mapping.MapDocument(mxd_path)
df  = arcpy.mapping.ListDataFrames(mxd)[0]  # primer Data Frame

# Capa + plantilla de labels
lyr = get_or_add_layer(df, shp_path, layer_name_in_toc)
tmpl = arcpy.mapping.Layer(lyr_template)
arcpy.mapping.UpdateLayer(df, lyr, tmpl, True)
try: lyr.showLabels = True
except: pass

where_base = base_where(expedientes)

# Bucle por rangos: centrado + escala ajustada a valores cerrados
for minimo, maximo, alias in rangos:
    where_rango = "{wb} AND {campo} > {minv} AND {campo} <= {maxv}".format(
        wb=where_base, campo=campo_rangos, minv=minimo, maxv=maximo
    )

    # 1) Filtro por rango
    lyr.definitionQuery = where_rango

    # 2) Asegurar estilo de etiquetas de la plantilla
    arcpy.mapping.UpdateLayer(df, lyr, tmpl, True)
    try: lyr.showLabels = True
    except: pass

    # 3) Calcular extensión del subconjunto
    tmp_view = "view_{}".format(alias)
    arcpy.MakeFeatureLayer_management(lyr.dataSource, tmp_view, where_rango)
    recs = int(arcpy.GetCount_management(tmp_view).getOutput(0))

    if recs > 0:
        desc = arcpy.Describe(tmp_view)
        if hasattr(desc, "extent") and desc.extent:
            ext = extent_expand(desc.extent, margen_ext)
            # Coloca la extensión y deja que ArcMap calcule una escala aproximada
            df.extent = ext
            escala_auto = df.scale
            # Ajusta la escala al valor "cerrado" inmediato superior
            escala_snapped = snap_escala(int(round(escala_auto)), escalas_cerradas)
            df.scale = escala_snapped
            # Recentrar (pan) con esa escala final para evitar cortes
            df.panToExtent(ext)

    arcpy.Delete_management(tmp_view)

    # 4) Exportar
    out_pdf = os.path.join(
        output_folder,
        u"Mapa_{0}_{1}regs_{2}.pdf".format(
            alias, recs, datetime.datetime.now().strftime("%Y%m%d")
        )
    )
    arcpy.mapping.ExportToPDF(mxd, out_pdf, df, resolution=dpi_export)
    print(u"Exportado: {} ({} registros)".format(out_pdf, recs))

# No guardamos el MXD (evita errores de permisos)
del mxd
print("Listo: PDFs con escala ajustada a valores 'cerrados' por rango.")

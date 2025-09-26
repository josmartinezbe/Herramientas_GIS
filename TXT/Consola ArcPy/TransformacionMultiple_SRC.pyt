# -*- coding: utf-8 -*-
import arcpy
import os

# ==========================
# PARÁMETROS (edítalos aquí)
# ==========================
in_path        = r"D:\drummond\SeleccionPredios_base.gdb\FDS_CRITERIOS_WGS84"  # FC o FDS de entrada
out_gdb        = r"D:\drummond\SeleccionPrediosCesar.gdb"                     # GDB de salida (debe existir)
out_fds_name   = "FDS_CRITERIOS"                                              # "" para salida en raíz de la GDB
sr_destino     = arcpy.SpatialReference(3116)                                 # 4326=WGS84, 9377=CTM12, etc.
transformacion = ""                                                           # Ej.: "MAGNA_Colombia_To_WGS_1984_7P"
prefijo        = ""                                                           # Ej.: "prj_"
sufijo         = ""                                                           # Ej.: "_wgs84"

arcpy.env.overwriteOutput = True

# ==========================
# FUNCIONES DE APOYO
# ==========================
def ensure_fds(gdb, fds_name, spatial_ref):
    """Crea el Feature Dataset si no existe. Devuelve ruta completa o la GDB si fds_name==""."""
    if not fds_name:
        return gdb
    fds_path = os.path.join(gdb, fds_name)
    if not arcpy.Exists(fds_path):
        arcpy.CreateFeatureDataset_management(gdb, fds_name, spatial_ref)
    else:
        # Validar que SR coincida con sr_destino
        try:
            if arcpy.Describe(fds_path).spatialReference.factoryCode != spatial_ref.factoryCode:
                raise RuntimeError(
                    "El SR del FDS de salida no coincide con el SR destino. "
                    "Crea un FDS con el SR correcto o cambia 'sr_destino'."
                )
        except Exception:
            pass
    return fds_path

def project_fc(in_fc, out_workspace, sr_dest, transf, prefix="", suffix=""):
    name = arcpy.Describe(in_fc).baseName
    out_fc = os.path.join(out_workspace, "{}{}{}".format(prefix, name, suffix))
    arcpy.Project_management(
        in_dataset=in_fc,
        out_dataset=out_fc,
        out_coor_system=sr_dest,
        transform_method=(transf or None)
    )
    arcpy.AddMessage("OK -> {}".format(out_fc))
    print("OK -> {}".format(out_fc))

# ==========================
# LÓGICA PRINCIPAL
# ==========================
if not arcpy.Exists(in_path):
    raise RuntimeError("No existe la entrada: {}".format(in_path))
if not arcpy.Exists(out_gdb):
    raise RuntimeError("No existe la GDB de salida: {}".format(out_gdb))

# Crear/validar FDS de salida (o usar raíz de la GDB)
out_workspace = ensure_fds(out_gdb, out_fds_name, sr_destino)

desc = arcpy.Describe(in_path)
in_dtype = desc.dataType  # "FeatureClass" o "FeatureDataset"

if in_dtype.lower() == "featureclass":
    # Caso 1: una sola capa
    project_fc(in_path, out_workspace, sr_destino, transformacion, prefijo, sufijo)

elif in_dtype.lower() == "featuredataset":
    # Caso 2: múltiples capas dentro de un FDS
    arcpy.env.workspace = in_path
    fcs = arcpy.ListFeatureClasses()
    if not fcs:
        raise RuntimeError("No se encontraron feature classes dentro de: {}".format(in_path))
    for fc in fcs:
        in_fc = os.path.join(in_path, fc)
        project_fc(in_fc, out_workspace, sr_destino, transformacion, prefijo, sufijo)
else:
    raise RuntimeError("La ruta de entrada debe ser un Feature Class o un Feature Dataset.")

# Save the geometry (triangles, verticies) of a FESOM grid to a gdal dataset
# Author: R. Rietbroek
# Date 17 May 2019
# Currently this save the surface nodes only
# Improvements are possible
# * Optionally store the 3D surfaces
# * add info on e.g. bathymetry  the nodes
# * 
import osgeo,ogr as ogr
import osgeo.osr as osr

def fesom2gdal(mesh,outputname,gdaldriver='GPKG'):
    """Export a FESOM surface mesh to a GIS shapefile
    input: mesh a FESOM mesh loaded with fesom_mesh(meshpath, get3d=True)
        outputname: the name of the output dataset 
        gdaldriver: the driver to use to write the output (defaults to geopackage, but could be anything the gdal library supports including POSTGIS)
    returns: nothing"""
    driver = ogr.GetDriverByName(gdaldriver)

    data_source = driver.CreateDataSource(outputname)

    # create the spatial reference, WGS84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)


    # create the layer containing the vertices
    vertlayer = data_source.CreateLayer("vert", srs, ogr.wkbPoint)

    field_onbound = ogr.FieldDefn("onboundary", ogr.OFTInteger)
    vertlayer.CreateField(field_onbound)
    
    field_topo = ogr.FieldDefn("topo", ogr.OFTReal)
    vertlayer.CreateField(field_topo)

    #also store the indices of the 3Delems for the corresponding z-levels (to look up 3d indices)
    # NOTE: ESRI SHapefiles do not support lists as attributes!! so this will not be registered when
    field_zindx=ogr.FieldDefn("nodeid",ogr.OFTIntegerList)
    vertlayer.CreateField(field_zindx)
    
    # add vertices 
    for id,(lon,lat,onbnd) in enumerate(zip(mesh.x2,mesh.y2,mesh.ind2d)):
        feature = ogr.Feature(vertlayer.GetLayerDefn())
        # Note: we need a conversion to int so the value get's accepted by the gdal library
        feature.SetField("onboundary", int(onbnd))

        feature.SetField('topo',mesh.topo[id])
        # note: we subtract 1 here to be consistent with the zero-indexing used in nodeid
        idxfield=feature.GetFieldIndex("nodeid")
        feature.SetFieldIntegerList(idxfield,[int(x-1) for x in mesh.n32[id,:] if x >=0])
        #create a point geometry
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(lon,lat)
        feature.SetGeometry(point)
        vertlayer.CreateFeature(feature)
        # Dereference the feature in order to appropriately call the destructor
        feature = None
    
    # create the layer containing the triangles
    # NOTE: the dedicated triangle type is not visible in Qgis
    # surftype=ogr.wkbTriangle
    surftype=ogr.wkbPolygon


    tinlayer = data_source.CreateLayer("tin", srs, surftype)

    # Add the fields we're interested in (nodeid's of 
    field_nodeid = ogr.FieldDefn("nodeid", ogr.OFTIntegerList)
    tinlayer.CreateField(field_nodeid)
    field_area = ogr.FieldDefn("area", ogr.OFTReal)
    tinlayer.CreateField(field_area)
    
    field_topo = ogr.FieldDefn("topo", ogr.OFTReal)
    tinlayer.CreateField(field_topo)

    # exclude cyclic elements
    elem2=mesh.elem[mesh.no_cyclic_elem,:]
    #loop over triangular elements
    for i1,i2,i3 in elem2:
        feature = ogr.Feature(tinlayer.GetLayerDefn())
        
        ring=ogr.Geometry(ogr.wkbLinearRing)
        tri=ogr.Geometry(surftype)
        
        ring.AddPoint(mesh.x2[i1],mesh.y2[i1])
        ring.AddPoint(mesh.x2[i2],mesh.y2[i2])
        ring.AddPoint(mesh.x2[i3],mesh.y2[i3])
        tri.AddGeometry(ring)
        
        idxfield=feature.GetFieldIndex("nodeid")
        feature.SetFieldIntegerList(idxfield,[int(i1),int(i2),int(i3)])
        # TODO convert to squared km (which projection is used for FESOM??)
        feature.SetField("area", tri.Area())
        
        #currently just set topo to the mean of the topo of the vertices
        feature.SetField("topo", (mesh.topo[i1]+mesh.topo[i2]+mesh.topo[i3])/3.0)
        
        feature.SetGeometry(tri)
        tinlayer.CreateFeature(feature)
        feature=None


    # Save and close the data source
    data_source = None


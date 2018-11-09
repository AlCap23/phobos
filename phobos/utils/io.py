#!/usr/bin/python
# coding=utf-8

# -------------------------------------------------------------------------------
# This file is part of Phobos, a Blender Add-On to edit robot models.
# Copyright (C) 2018 University of Bremen & DFKI GmbH Robotics Innovation Center

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------------

import shutil
import sys
import os
import bpy
import numpy

from phobos import defs
from phobos import display
from phobos.phoboslog import log

from phobos.io.entities import entity_types
from phobos.io.meshes import mesh_types
from phobos.io.scenes import scene_types

from phobos.utils import selection as sUtils
from phobos.utils import naming as nUtils
from phobos.utils import blender as bUtils



indent = '  '
xmlHeader = '<?xml version="1.0"?>\n<!-- created with Phobos ' + defs.version + ' -->\n'


def xmlline(ind, tag, names, values):
    """Generates an xml line with specified values.
    To use this function you need to know the indentation level you need for this line.
    Make sure the names and values list have the correct order.

    Args:
      ind(int >= 0): Indentation level
      tag(String): xml element tag
      names(list (same order as for values): Names of xml element's attributes
      values(list (same order as for names): Values of xml element's attributes

    Returns:
      : String -- Generated xml line.

    """
    line = [indent * max(0, ind) + '<' + tag]
    for i in range(len(names)):
        line.append(' ' + names[i] + '="' + str(values[i]) + '"')
    line.append('/>\n')
    return ''.join(line)


def l2str(items, start=0, end=-1):
    """Generates string from (part of) a list.

    Args:
      items(list): List from which the string is derived (elements need to implement str())
      start(int, optional): Inclusive start index for iteration (Default value = 0)
      end(int, optional): Exclusive end index for iteration (Default value = -1)

    Returns:
      : str - Generated string.

    """
    start = max(start, 0)
    end = end if end >= 0 else len(items)
    return ' '.join([str(i) for i in items[start:end]])


def securepath(path):
    """Checks whether directories of a path exist and generates them if necessary.

    Args:
      path(str): The path to be secured (directories only)

    Returns:
      : String -- secured path as absolute path, None on error

    """
    path = os.path.abspath(path)
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except NotADirectoryError:
            log(path + " is not a valid directory", "ERROR")
            return None
    return path


def getExpSettings():
    """Returns Phobos' export settings as displayed in the GUI"""
    return bpy.context.scene.phobosexportsettings


def getExportModels():
    """Returns a list of objects representing a model (root) in the Blender scene"""
    if getExpSettings().selectedOnly:
        roots = [root for root in sUtils.getRoots() if root.select]
    else:
        roots = sUtils.getRoots()
    return list(roots)


def getEntityRoots():
    """Returns all objects for which entity properties are specified.

    Args:

    Returns:
      list: roots of well-defined entities in the scene

    """
    roots = [
        obj
        for obj in bpy.context.scene.objects
        if sUtils.isEntity(obj) and (not getExpSettings().selectedOnly or obj.select)
    ]
    return roots


def getModelListForEnumProp(self, context):
    """Returns list of all exportable models in the scene formatted for an EnumProperty

    Args:
      context: 

    Returns:

    """
    rootnames = set([nUtils.getModelName(r) for r in getExportModels()])
    return sorted([(r,) * 3 for r in rootnames])


def getDictFromYamlDefs(phobostype, defname, name):
    """Returns a phobos representation of the object specified by the definition parameters.

    Args:
      phobostype(str): phobostype of the definition
      defname(str): name of the individual definition
      name(str): name for the new object

    Returns:
      : dict -- phobos representation from the definition

    """
    if 'material' in defs.def_settings[phobostype + 's'][defname]:
        material = defs.def_settings[phobostype + 's'][defname]['material']
    else:
        material = None

    # separate properties and annotations from each other
    props = {
        key: value
        for key, value in defs.definitions[phobostype + 's'][defname].items()
        if not isinstance(value, dict)
    }
    annots = {
        key: value
        for key, value in defs.definitions[phobostype + 's'][defname].items()
        if isinstance(value, dict)
    }

    # create a dictionary holding the properties
    phobos_dict = {
        'name': name,
        'defname': defname,
        'category': defs.def_settings[phobostype + 's'][defname]['categories'],
        'material': material,
        'type': defs.def_settings[phobostype + 's'][defname]['type'],
        'props': props,
        'annotations': annots,
    }

    # add the general settings for this object
    general_settings = defs.def_settings[phobostype + 's'][defname]
    for gensetting in general_settings.keys():
        phobos_dict[gensetting] = general_settings[gensetting]

    return phobos_dict


def getOutputMeshtype():
    """Returns the mesh type to be used in exported files as specified in the GUI"""
    return str(getExpSettings().outputMeshtype)


def getOutputMeshpath(path, meshtype=None):
    """Returns the folder path for mesh file export as specified in the GUI.
    
    Phobos by default creates a directory 'meshes' in the export path and subsequently creates
    sub-directories of "export/path/meshes" for every format, e.g. resulting in "export/path/mesh/obj"
    for .obj file export.

    Args:
      path(str): export path root (set in the GUI)
      meshtype(str, optional): a valid mesh type, otherwise the type set in the GUI is used (Default value = None)

    Returns:
      string: output path for meshes

    """
    return os.path.join(path, 'meshes', meshtype if meshtype else getOutputMeshtype())


def getEntityTypesForExport():
    """Returns list of entity types available for export"""
    return [
        typename
        for typename in sorted(entity_types)
        if 'export' in entity_types[typename] and 'extensions' in entity_types[typename]
    ]


def getEntityTypesForImport():
    """Returns list of entity types available for import"""
    return [
        typename
        for typename in sorted(entity_types)
        if 'import' in entity_types[typename] and 'extensions' in entity_types[typename]
    ]


def getSceneTypesForExport():
    """Returns list of scene types available for export"""
    return [
        typename
        for typename in sorted(scene_types)
        if 'export' in scene_types[typename] and 'extensions' in scene_types[typename]
    ]


def getSceneTypesForImport():
    """Returns list of scene types available for import"""
    return [
        typename
        for typename in sorted(scene_types)
        if 'import' in scene_types[typename] and 'extensions' in scene_types[typename]
    ]


def getMeshTypesForExport():
    """Returns list of mesh types available for export"""
    return [
        typename
        for typename in sorted(mesh_types)
        if 'export' in mesh_types[typename] and 'extensions' in mesh_types[typename]
    ]


def getMeshTypesForImport():
    """Returns list of mesh types available for import"""
    return [
        typename
        for typename in sorted(mesh_types)
        if 'import' in mesh_types[typename] and 'extensions' in mesh_types[typename]
    ]


def getExportPath():
    """Returns the root path of the export.
    
    Phobos by default creates directories in this path for every format that is used for export,
    i.e. a folder "export/path/urdf" will be created for URDF files.

    Args:

    Returns:

    """
    if os.path.isabs(getExpSettings().path):
        return getExpSettings().path
    else:
        return os.path.normpath(os.path.join(bpy.path.abspath('//'), getExpSettings().path))


def getAbsolutePath(path):
    """Returns an absolute path derived from the .blend file

    Args:
      path: 

    Returns:

    """
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(bpy.path.abspath('//'), path)


def importBlenderModel(filepath, namespace='', prefix=False):
    """Imports an existing Blender model into the current .blend scene

    Args:
      filepath(str): Path of the .blend file
      namespace: (Default value = '')
      prefix: (Default value = False)

    Returns:

    """
    if os.path.exists(filepath) and os.path.isfile(filepath) and filepath.endswith('.blend'):
        log("Importing Blender model" + filepath, "INFO")
        objects = []
        with bpy.data.libraries.load(filepath) as (data_from, data_to):
            for objname in data_from.objects:
                objects.append({'name': objname})
        bpy.ops.wm.append(directory=filepath + "/Object/", files=objects)
        imported_objects = bpy.context.selected_objects
        resources = [obj for obj in imported_objects if obj.name.startswith('resource::')]
        new_objects = [obj for obj in imported_objects if not obj.name.startswith('resource::')]
        if resources:
            if 'resources' not in bpy.data.scenes.keys():
                bpy.data.scenes.new('resources')
            sUtils.selectObjects(resources)
            bpy.ops.object.make_links_scene(scene='resources')
            bpy.ops.object.delete(use_global=False)
            sUtils.selectObjects(new_objects)
        bpy.ops.view3d.view_selected(use_all_regions=False)
        # allow the use of both prefixes and namespaces, thus truly merging
        # models or keeping them separate for export
        if namespace != '':
            if prefix:
                for obj in bpy.context.selected_objects:
                    # set prefix instead of namespace
                    obj.name = namespace + '__' + obj.name
                    # make sure no internal name-properties remain
                    for key in obj.keys():
                        try:
                            if obj[key].endswidth("/name"):
                                del obj[key]
                        except AttributeError:  # prevent exceptions from non-string properties
                            pass
            else:
                for obj in bpy.context.selected_objects:
                    nUtils.addNamespace(obj, namespace)
        submechanism_roots = [
            obj
            for obj in bpy.data.objects
            if obj.phobostype == 'link' and 'submechanism/spanningtree' in obj
        ]
        for root in submechanism_roots:
            partlist = [root] + root['submechanism/spanningtree']
            if 'submechanism/freeloader' in root:
                partlist += root['submechanism/freeloader']
            sUtils.selectObjects(partlist, active=0)
            bpy.ops.group.create(name='submechanism:' + root['submechanism/name'])
        return True
    else:
        return False


def importResources(restuple, filepath=None):
    """Accepts an iterable of iterables describing resource objects to import. For instance,
    reslist=(('joint', 'continuous'), ('interface', 'default', 'bidirectional'))
    would import the resource objects named 'joint_continuous' and
    'interface_default_bidirectional'.

    Args:
      restuple: iterable of iterables containing import objects
      filepath: path to file from which to load resource (Default value = None)
    Returns(tuple of bpy.types.Object): imported objects

    Returns:

    """
    currentscene = bpy.context.scene.name
    bUtils.switchToScene('resources')
    # avoid importing the same objects multiple times
    imported_objects = [
        nUtils.stripNamespaceFromName(obj.name) for obj in bpy.data.scenes['resources'].objects
    ]
    requested_objects = ['_'.join(resource) for resource in restuple]
    new_objects = [obj for obj in requested_objects if obj not in imported_objects]

    # if no filepath is provided, use the path from the preferences
    if not filepath:
        filepath = os.path.join(bUtils.getPhobosConfigPath(), 'resources', 'resources.blend')
        print(filepath)

    # import new objects from resources.blend
    if new_objects:
        with bpy.data.libraries.load(filepath) as (data_from, data_to):
            objects = [{'name': name} for name in new_objects if name in data_from.objects]
            if objects:
                bpy.ops.wm.append(directory=filepath + "/Object/", files=objects)
            else:
                log('Resource objects could not be imported.', 'ERROR')
                bUtils.switchToScene(currentscene)
                return None
    objects = bpy.context.selected_objects
    for obj in objects:
        nUtils.addNamespace(obj, 'resource')
    bUtils.switchToScene(currentscene)
    return objects


def getResource(specifiers):
    """Returns a resource object defined by an iterable of strings.

    Args:
      specifiers(iterable): strings specifying the resource

    Returns:
      : bpy.types.Object -- resource object (or None if it could not be imported)

    """
    log("Searching for resource object " + '_'.join(specifiers) + ".", 'DEBUG')

    resobjname = nUtils.addNamespaceToName('_'.join(specifiers), 'resource')

    if 'resources' not in bpy.data.scenes or resobjname not in bpy.data.scenes['resources'].objects:
        newobjects = importResources((specifiers,))
        if not newobjects:
            return None
        return newobjects[0]
    return bpy.data.scenes['resources'].objects[resobjname]


def copy_model(model):
    """Returns a recursive deep copy of a model dictionary.
    
    The deep copy recreates dictionaries and lists, while keeping Blender objects and everything
    else untouched.
    
    This function is required, as we can not use copy.deepcopy() due to the Blender objects in our
    Phobos representation.

    Args:
      model(dict): model dictionary to copy

    Returns:
      : dict -- deep copy of the model dictionary

    """
    if isinstance(model, dict):
        newmodel = {}
        for key, value in model.items():
            if isinstance(value, dict) or isinstance(value, list):
                newmodel[key] = copy_model(value)
            else:
                newmodel[key] = value
        return newmodel
    elif isinstance(model, list):
        newlist = []
        for value in model:
            if isinstance(value, bpy.types.Object):
                newlist.append(value)
            elif isinstance(value, dict) or isinstance(value, list):
                newlist.append(copy_model(value))
            else:
                newlist.append(value)
        return newlist
    raise TypeError(
        "Deep copy failed. Unsuspected element in the dictionary: {}".format(type(model))
    )


def exportModel(model, exportpath='.', entitytypes=None):
    """Exports model to a given path in the provided formats.

    Args:
      model(dict): dictionary of model to export
      exportpath(str, optional): path to export root (Default value = '.')
      entitytypes(list of str, optional): export types - model will be exported to all (Default value = None)

    Returns:

    """
    if not exportpath:
        exportpath = getExportPath()
    if not entitytypes:
        entitytypes = getEntityTypesForExport()

    # TODO: Move texture export to individual formats? This is practically SMURF
    # TODO: Also, this does not properly take care of textures embedded in a .blend file
    # export textures
    if getExpSettings().exportTextures:
        path = os.path
        for materialname in model['materials']:
            mat = model['materials'][materialname]
            for texturetype in ['diffuseTexture', 'normalTexture', 'displacementTexture']:
                # skip materials without texture
                if texturetype not in mat:
                    continue

                sourcepath = path.join(path.expanduser(bpy.path.abspath('//')), mat[texturetype])
                if path.isfile(sourcepath):
                    texture_path = securepath(path.join(exportpath, 'textures'))
                    log(
                        "Exporting texture {} of material {} to {}.".format(
                            texturetype, mat[texturetype], texture_path
                        ),
                        'INFO',
                    )
                    try:
                        shutil.copy(
                            sourcepath, path.join(texture_path, path.basename(mat[texturetype]))
                        )
                    except shutil.SameFileError:
                        log(
                            "{} of material {} already in place.".format(texturetype, materialname),
                            'WARNING',
                        )

                    # update the texture path in the model
                    mat[texturetype] = 'textures/' + path.basename(mat[texturetype])

    # export model in selected formats
    for entitytype in entitytypes:
        typename = "export_entity_" + entitytype
        # check if format exists and should be exported
        if not getattr(bpy.context.scene, typename, False):
            continue
        # format exists and is exported:
        model_path = os.path.join(exportpath, entitytype)
        securepath(model_path)

        # export model using entity export function
        log("Export model '" + model['name'] + "' as " + entitytype + " to " + model_path, "DEBUG")

        # pass a model copy to the entity export, as these might alter the dictionary
        newmodel = copy_model(model)
        entity_types[entitytype]['export'](newmodel, model_path)

    # export meshes in selected formats
    i = 1
    mt = len([m for m in mesh_types if getattr(bpy.context.scene, "export_mesh_" + m, False)])
    mc = len(model['meshes'])
    n = mt * mc
    for meshtype in mesh_types:
        mesh_path = getOutputMeshpath(exportpath, meshtype)
        try:
            if getattr(bpy.context.scene, "export_mesh_" + meshtype, False):
                securepath(mesh_path)
                for meshname in model['meshes']:
                    mesh_types[meshtype]['export'](model['meshes'][meshname], mesh_path)
                    display.setProgress(i / n, 'Exporting ' + meshname + '.' + meshtype + '...')
                    i += 1
        except KeyError as e:
            log("Error exporting mesh {0} as {1}: {2}".format(meshname, meshtype, str(e)), "ERROR")
    display.setProgress(0)


def exportScene(
    scenedict, exportpath='.', scenetypes=None, export_entity_models=False, entitytypes=None
):
    """Exports provided scene to provided path

    Args:
      scenedict(dict): dictionary of scene
      exportpath(str, optional): path to scene export folder (Default value = '.')
      scenetypes(list of str, optional): export types for scene - scene will be exported to all (Default value = None)
      export_entity_models(bool, optional): whether to export entities additionally (Default value = False)
      entitytypes(list of str, optional): types to export entities in in case they are exported (Default value = None)

    Returns:

    """
    if not exportpath:
        exportpath = getExportPath()
    if not scenetypes:
        scenetypes = getSceneTypesForExport()
    if export_entity_models:
        for entity in scenedict['entities']:
            exportModel(entity, exportpath, entitytypes)
    for scenetype in scenetypes:
        gui_typename = "export_scene_" + scenetype
        # check if format exists and should be exported
        if getattr(bpy.context.scene, gui_typename):
            scene_types[scenetype]['export'](
                scenedict['entities'], os.path.join(exportpath, scenedict['name'])
            )


def import_csv(filepath, skip_header = 0, delimiter = ","):
    """Imports a csv file and returns a dictionary.
    """

    # Import the csv
    data_dict = {}
    try:
        data = numpy.genfromtxt(filepath, delimiter = delimiter, skip_header = skip_header, names = True)
        names = data.dtype.names
        for name in names:
            data_dict.update({name : list(data[name])})
    except KeyError as err:
        log("Error : {} not a valid csv file.".format(filepath), level = "ERROR")
        return {}

    return data_dict

def create_Reachability(obj, column_dict, name, size, shape):
    """Creates a reachability map from a given csv data object.
    """

    from phobos.utils.editing import parentObjectsTo

    # Dict should have structure like:
    # x : Columnname, y: columnname ...
    # m : Measurement, e.g. Dexterity measure

    # Gather the coordinates
    vertices = []
    for i,coordinates in enumerate(['x', 'y', 'z']):
        try:
            vertices.append(obj[column_dict[coordinates]])
        except KeyError as err:
            log("Error : {} not a valid key inside the supplied data.".format(coordinates), level = "ERROR")

    vertices = numpy.array(vertices)
    coordinates = []

    for i in range(vertices.shape[1]):
        coordinates.append(tuple(vertices[:, i]))
    

    # Start the visualization

    # Create a mesh
    mesh = bpy.data.meshes.new(name)
    mesh.update()
    mesh.validate()

    mesh.from_pydata(coordinates, [], [])
    
    # Add a mesh obj
    scene = bpy.context.scene
    mesh_obj = bpy.data.objects.new(name, mesh)
    mesh_obj['phobostype'] = 'data'
    for i in range(20):
        mesh_obj.layers[i] = ( i == defs.layerTypes["data"])
    scene.objects.link(mesh_obj)
    sUtils.selectObjects([mesh_obj])

    # Add the visualization obj
    #bpy.ops.mesh.primitive_uv_sphere_add(size = 0.0025)
    mesh_name = nUtils.getUniqueName(name + '_Primitve', bpy.data.objects)
    if shape == 'sphere':
        viz_mesh = bUtils.createPrimitive(
            mesh_name,
            'sphere',
            size,
            defs.layerTypes["data"],
            phobostype = 'data',
        )
    else:
        viz_mesh = bUtils.createPrimitive(
            mesh_name,
            'box',
            (size, ) * 3,
            defs.layerTypes["data"],
            phobostype = 'data',
        )
    

    # Add a particle system
    material_name = "Visualization_Mat"
    material = bpy.data.materials.new(name = material_name)
    material.diffuse_color = (1, 0, 0) 
    material.use_transparency = True
    viz_mesh.data.materials.append(material)

    # Add the material in cylces
    bpy.context.scene.render.engine = 'CYCLES'
    material.use_nodes = True
    node_tree = material.node_tree

    if 'Material Output' in node_tree.nodes:
        material_output_node = node_tree.nodes['Material Output']
    else:
        material_output_node = node_tree.nodes.new('ShaderNodeOutputMaterial')
    
    if 'Diffuse BSDF' in node_tree.nodes:
        diffuse_node = node_tree.nodes['Diffuse BSDF']
    else:
        diffuse_node = node_tree.nodes.new("ShaderNodeBsdfDiffuse")
    node_tree.links.new(diffuse_node.outputs['BSDF'], material_output_node.inputs['Surface'])
    
    if 'Image Texture' in node_tree.nodes:
        image_texture_node = node_tree.nodes['Image Texture']
    else:
        image_texture_node = node_tree.nodes.new("ShaderNodeTexImage")
    node_tree.links.new(image_texture_node.outputs['Color'], diffuse_node.inputs['Color'])
    
    # To view the texture we set the height of the texture to vis_image_height 
    vis_image_height = 1
    image = bpy.data.images.new('ParticleColor', len(coordinates), vis_image_height)
    local_pixels = list(image.pixels[:])
    num_points = vertices.shape[1]
    
    
    # Normalize the measurement
    if column_dict['w']:
        try:
            measurement = numpy.array(obj[column_dict['w']])
            min_val = numpy.max(measurement)
            max_val = numpy.min(measurement)
            normalized = (measurement - min_val)/(max_val - min_val)
        except:
            log("No measurement data found!", level = "ERROR")
            pass
        

        for j in range(vis_image_height):
            for point_index, point in enumerate(coordinates):
                column_offset = point_index * 4
                row_offset = j * 4 * num_points
                local_pixels[row_offset + column_offset] =  0 if normalized[point_index] <= 0.5 else 2.*(normalized[point_index] - 0.5)
                local_pixels[row_offset + column_offset + 1] =  2.*normalized[point_index] if normalized[point_index] <= 0.5 else 1.-2.*(normalized[point_index]-0.5) 
                local_pixels[row_offset + column_offset + 2] = 2.*(0.5-normalized[point_index]) if normalized[point_index] <= 0.5 else 0
        image.pixels = local_pixels[:]
    
    image_texture_node.image = image
    particle_info_node = node_tree.nodes.new('ShaderNodeParticleInfo')
    divide_node = node_tree.nodes.new('ShaderNodeMath')
    divide_node.operation = 'DIVIDE'
    node_tree.links.new(particle_info_node.outputs['Index'], divide_node.inputs[0])
    divide_node.inputs[1].default_value = num_points
    shader_node_combine = node_tree.nodes.new('ShaderNodeCombineXYZ')
    node_tree.links.new(divide_node.outputs['Value'], shader_node_combine.inputs['X'])
    node_tree.links.new(shader_node_combine.outputs['Vector'], image_texture_node.inputs['Vector'])

    if len(mesh_obj.particle_systems) == 0:
        mesh_obj.modifiers.new("particle sys", type='PARTICLE_SYSTEM')
        particle_sys = mesh_obj.particle_systems[0]
        settings = particle_sys.settings
        settings.type = 'HAIR'
        settings.use_advanced_hair = True
        settings.emit_from = 'VERT'
        settings.count = len(coordinates)
        settings.hair_length = 100           # This must not be 0
        settings.use_emit_random = False
        settings.render_type = 'OBJECT'
        settings.dupli_object = viz_mesh

    parentObjectsTo([viz_mesh], mesh_obj, clear = True)
    return mesh_obj

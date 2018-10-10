#!/usr/bin/python
# coding=utf-8

import os
import shutil

import bpy
import mathutils

import phobos.defs as defs
from phobos.utils.io import xmlHeader
from phobos.utils.io import indent as phobosindentation
from phobos.utils.io import l2str as list_to_string
from phobos.utils.io import getExpSettings
import phobos.utils.general as gUtils
from phobos.utils.blender import getPhobosPreferences
import phobos.utils.selection as sUtils
import phobos.model.models as models
from phobos.phoboslog import log

class luaTagger(object):
    """A simple class to create a model in lua format.

    Args:

    Returns:
    """

    def __init__(self, indentation = 0, initial = 0):
        """ Creates a new lua tagger. The indentation is used for nesting of vectors.

            :param indentation: level of hierarchy
            :type indentation: int
            :param initial: initial level of hierarchy
            :type initial: int
        """

        self.indentation = indentation
        self.initial = initial
        self.workingTags = []
        self.output = ""
    
    def ind(self, add = 0):
        """Helper function to return the current indentation depending on the
        hierarchy.
        
        :return: str -- the current indentation (e.g. "  ").

        Args:

        Returns:

        """
        return "" + (self.indentation + add ) * "    "

    def ascend(self):
        """Move up one hierarchical layer by finishing the current tag and
        removing one indentation.
        
        :exception IndentationError -- trying to move above root hierarchical
        layer

        Args:

        Returns:

        """
        if self.indentation > self.initial:
            lasttag = self.workingTags.pop(-1)
            self.indentation -= 1
            self.output += "{\n" + lasttag + "\n}"
        else:
            IndentationError()

    def descend(self, tag = None, attribs=None):
        """Move down one hierarchical layer using the new tag.
        Optional in-line attributes can be provided in the
        dictionary attribs (e.g. {'name': 'foo'}.

        Args:
          tag(str): tag to descend with
          attribs: (Default value = None)

        Returns:
          None: None

        """

        # Store special tags
        if tag is not None:
            self.workingTags.append(str(tag))
        
        # create parameter strings by unpacking dictionary
        if attribs:
            parameters = [key + '= {\n' + self.ind(1) + str(attribs[key]) + ' \n}, ' for key in attribs.keys()]
            # remove trailing whitespace
            parameters[-1] = parameters[-1][:-1]
        else:
            parameters = {}

        line = self.ind() + tag + (" = { \n" if len(parameters) > 0 else "")

        # add optional parameters
        if len(parameters) > 0:
            for param in parameters:
                line += param
            line += "\n"+ self.ind() + "},"

        # finish line and descend one layer
        self.output += line + "\n"
        self.indentation += 1


    def write(self, text):
        """Write a custom line to the output. Use to create the header or comments.

        Args:
          text(str): The line to write (line break has to be included)

        Returns:

        """
        self.output += text

    def valueVarTable(self, name, values):
        """Adds an attribute to the current element. The tag of the attribute
        is wrapped around its value.

        Args:
        name(str) : ID of the table
          tags(list): The tags of the attributes.
          value(list): The list of values of the attribute

        Returns:

        """
        value_string = self.ind() + str(name) + ' = {\n'
        
        for tag, value in zip(tags,values):
                value_string += self.ind(1) + '\'' + str(value) + '\', \n'
        value_string += self.ind() + '}'

        self.output += value_string

    def valueTable(self, name, tags, values):
        """Adds an attribute to the current element. The tag of the attribute
        is wrapped around its value.

        Args:
        name(str) : ID of the table
          tags(list): The tags of the attributes.
          value(list): The list of values of the attribute

        Returns:

        """
        value_string = self.ind() + str(name) + ' = {\n'
        
        for tag, value in zip(tags,values):
                value_string += self.ind(1) + str(tag) +' = ' + str(value) + ', \n'
        value_string += self.ind() + '}'

        self.output+= value_string

    def value(self, tag, value, table = False):
        """Adds an attribute to the current element. The tag of the attribute
        is wrapped around its value.

        Args:
          tag(str): The tag of the attribute.
          value(str (will be casted anyway): The value of the attribute

        Returns:

        """
        if isinstance(value, list) and value:
            str_value = "{" 
            if table:
                str_value += "\n{ "
            if len(value) > 1 :
                for val in value[:-1]:
                    str_value += " {0},".format(val)
                str_value += " {} }}".format(value[-1])
            else:
                str_value += " {} }}".format(value[0])
            
            if table:
                str_value += "\n}"
        elif isinstance(value, list) and not value:
            str_value = ""
            if table:
                str_value += "{\n"
            str_value += "{ }"
            if table:
                str_value += " }"
        else:
            str_value = ""
            if table:
                str_value += "{ "
            str_value += str(value)
            if table:
                str_value += "\n}"

        self.output += self.ind() + str(tag) + '= ' + str(str_value) + ',\n'

    def valueVar(self, value):
        """Adds an attribute to the current element. The tag of the attribute
        is wrapped around its value.
        Args:
          tag(str): The tag of the attribute.
          value(str (will be casted anyway): The value of the attribute
        Returns:
        """
        self.output += self.ind() + '= \'' + str(value) + '\', \n'

    def get_indent(self):
        """TODO Missing documentation"""
        return self.indentation

    def get_output(self):
        """Completes all trailing tags until at initial indentation and
        returns the output as string.
        
        :return: str -- the finished xml string.

        Args:

        Returns:

        """
        # ascend to base layer
        while self.indentation > self.initial:
            self.ascend()

        return self.output


def exportRBDLInertial(indentation, name, inertialdata):
    """Creates an inertia table
    
    Args:
     indentation(int) : indentation at current level
     name(str): Name of the inertia
     mass(float): Mass value in kg
     com(list): Coordinates of the center of mass in parent frame
     inertia(list): Inertia of the body given at the center of mass

    Returns:
     string representation of the inertial
    """

    tagger = luaTagger(initial = indentation)
    
    if inertialdata:
        print(inertialdata.keys())
        
        mass = inertialdata['mass']
        
        # We always export at the link
        com = [0., 0., 0.]

        # Preprocess the inertia
        i_xx = inertialdata['inertia'][0] 
        i_xy = inertialdata['inertia'][1]
        i_xz = inertialdata['inertia'][2]
        i_yy = inertialdata['inertia'][3]
        i_yz = inertialdata['inertia'][4]
        i_zz = inertialdata['inertia'][5]
    else:
        mass = 1e-3
        com = [0., 0., 0.]

        i_xx = 1e-3
        i_xy = 0.0
        i_xz = 0.0
        i_yy = 1e-3
        i_yz = 0.0
        i_zz = 1e-3

    inertia_table = '{ \n'
    inertia_table += '  {{ {0}, {1}, {2} }}, \n'.format(i_xx, i_xy, i_xz)
    inertia_table += '  {{ {0}, {1}, {2} }}, \n'.format(i_xy, i_yy, i_yz)
    inertia_table += '  {{ {0}, {1}, {2} }}, \n'.format(i_xz, i_yz, i_zz)
    inertia_table += '}'

    # Preprocess the com
    com_table = '{{ {0}, {1}, {2} }}'.format(com[0], com[1], com[2])

    tagger.valueTable(name, ['mass', 'com', 'inertia'], [mass, com_table, inertia_table])

    return "".join(tagger.get_output())

def exportRBDLPose(indentation, name,  pose = None, poseobject = None):
    """
    """

    tagger = luaTagger(initial = indentation)

    if not poseobject:
        if pose:
            matrix = mathutils.Matrix(pose['rawmatrix'])
        else:
            posedata = {'r' : '{ 0., 0., 0.}', 'E' : '{ {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.} }' }
    else:
        # We take the pose w.r.t. the parent
        matrix = poseobject.matrix_local
    
    translation = matrix.to_translation()
    rotation = matrix.to_3x3()
    
    posedata = {'r' : list(translation), 
                'E' : [list(vector) for vector in list(rotation)]
                }
    
    posedata = gUtils.roundFloatsInDict(posedata, getExpSettings().decimalPlaces)

    translation = '{{ {0}, {1}, {2} }}'.format(posedata['r'][0], posedata['r'][1], posedata['r'][2])
    
    rotation_table = '{ \n'
    for row in posedata['E']:
        rotation_table +='  {'
        for col in row:
            rotation_table += ' {},'.format(col)
        rotation_table += ' }, \n'
    rotation_table += '} \n'

    tagger.valueTable(name, ['r', 'E'], [translation, rotation_table])

    return "".join(tagger.get_output())

def exportRBDLJoint(indentation, joint):
    """
    """

    name = joint['name']
    jointtype = joint['type']
    
    
    # Define a joint mapping of primitves
    tagger = luaTagger(initial = indentation)

    
    if jointtype == 'fixed':
        axis = []
    
    elif jointtype == 'prismatic':
        axis = [0., 0., 0., ] + list(joint['axis'])
        
    elif jointtype == 'revolute':
        axis =  list(joint['axis']) + [0., 0., 0., ] 
        
    else:
        log("No valid joint type found for {0}".format(name), 'ERROR')
        return ""
    

    log("Parsed joint {0} of type {1}".format(name, jointtype), 'INFO')   

    tagger.value(name, axis, table = True) 

    return "".join(tagger.get_output())
    

def exportRBDLFrame(indentation, name, parent, joint = None, joint_frame = None, bodies = None):
    """Creates a frame from the given parameter.

    Args:
    indentation(int) : indentation at current level
      name(str): Name of the frame
      parent(str): Name of the parent frame
      joint(list) : Specifies the type of joint
      joint_frame(list) : Specifies the origin of the joint in the frame of the parent
      body(list): Specification of the dynamical parameters of the body

    Returns:
      string representation of the frame in lua
    """

    tagger = luaTagger(initial = indentation)

    tagger.value('name', name)
    tagger.value('parent', parent)
    if joint:
        tagger.value('joint', joint)
        if joint_frame:
            tagger.value('joint_frame', joint_frame)
    if bodies:
        tagger.valueVarTable('bodies', bodies)
    
    return "".join(tagger.get_output())



    


def exportLUA(model, filepath):
    """
    """
    log("Export LUA to {}.".format(filepath), 'INFO')
    filename = os.path.join(filepath, model['name'] + '.lua')

    annotationdict = models.gatherAnnotations(model)

    log("Exporting '{0}'...".format(model['name']), "DEBUG")
    print('Model links implemented.')
    print('Model joints implemented.')
    print('Model visuals implemented.')

    lua = luaTagger()

    # Init the stuff
    joints = luaTagger()
    inertials = luaTagger()
    frames = luaTagger()
    pose = luaTagger()

    log("Parsing joints...", 'INFO')
    for jointkey in model['joints'].keys():
        joint = model['joints'][jointkey]
        joints.write(exportRBDLJoint(lua.get_indent(1), joint))
    
    log("Parsing links...", 'INFO')
    for linkkey in model['links'].keys():
        link = model['links'][linkkey]
        # Write the pose
        pose.write(exportRBDLPose(joints.get_indent(1), link['name'] + '_pose', pose = link['pose'], poseobject = None))
        # Write the inertial
        #print(link.keys())
        inertials.write(exportRBDLInertial(joints.get_indent(1), link['name'] + '_inertia', link['inertial']) )
    
    lua.value("joints", joints.get_output(), table=True)
    lua.value("inertials", inertials.get_output(), table=True)
    lua.value("poses", pose.get_output(), table=True)
    #lua.write(inertials.get_output()[0])
    #print(lua.get_output())
    
    outputtext = lua.get_output()

    log("Writing output...", 'INFO')
    with open(filename, 'w') as outputfile:
            outputfile.writelines(outputtext)
    
    #print(lua.get_output())

# registering export functions of types with Phobos
entity_type_dict = {'lua': {'export': exportLUA, 'extensions': ('lua')}}

    





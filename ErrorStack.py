# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 03:02:46 2021

@author: Filip
"""
import numpy as np
import hyperspy.api as hs
from specerr import SpecErrAbs2D
from specMapDiff import *



class ErrorStack():
    '''
     Creates a stack och ErrorMaps based on the number of components that can be matched between the original sought after component and NMF/BSS output
     From the input of factors and loadings components are reconstructed and compared with the original simulated component which are given in originalComponents. 
     With the parameter order the order of the original components can be given this will allow for the creation of matching metadata after the matching of components is done
     
     Which elements and corresponding x-ray lines to filter the relative error maps with is also taken as an optional argument with Ag and Cu being default 
     
     The stack of errormaps can be found in .stack
     .stack contains a list of all ErrorMaps, one for each component
     '''
     
    def __init__(self, facs,loads,originalComponents,order = None, elements = ['Ag', 'Cu'], xlines = ['Ag_La','Cu_Ka']):
        self.facs = facs
        self.loads = loads
        self.originalComponents = originalComponents
        self.elements = elements
        self.xlines = xlines
        self.order = order
        self.createErrorMaps()
                 
    def __repr__(self):
        return "<ErrorStack, Components: %d of types: %s with elements %s>" %(self.size, self.order, self.elements)
    
    def createErrorMaps(self):
        '''
        Optional parameter: order = list with names of orignalComponents
        Takes factor and loadings from NMF/BSS and compares them to the orignial simulated components 
        After matches are done based on AbsErr2D errormaps are created and returned
        
        '''
        
        self.reconstructed = cLoadsFacs(self.loads,self.facs)
        rcDim = len(self.reconstructed)
        ocDim = len(self.originalComponents)
        
        NavDiff = rcDim-ocDim
        
        #Warning that data is lost is dimensions doesn't match betweent the  reconstructed and original
        if (NavDiff > 0 ):
            print(str(NavDiff)+' more reconstruced dimention than original these will be discarded')
        elif (NavDiff < 0):
            print('One or more original components that reconstructions, these will not be included')
        self.size = ocDim
        self.stack = []
        
        
        #For each reconstructed component the absolute error is compared to the original components 
        #This is done using SpecErrAbs2D and the lowest error is saved as a match
        #Möjlig TODO: Lägga till någon typ av cutoff om man är ute efter fler komponenter ex. vid över 50% diff ta det inte som match
        
        for rc in self.reconstructed:
            bestIndex = 0
            bestValue = SpecErrAbs2D(rc,self.originalComponents[0])
            for j in range(1,ocDim):
                newValue = SpecErrAbs2D(rc,self.originalComponents[j])
                if (newValue<bestValue):
                    bestValue = newValue
                    bestIndex = j
            self.stack.append(ErrorMap(rc,self.originalComponents[bestIndex],elements = self.elements, xlines = self.xlines, componentType = self.order[bestIndex]))
                
    def plot(self, component = 'all', errorType = 'all', elements = 'all'):
        '''
        Väljer vilka componenter som ska plottas, lite lurigare om man kan ha 
        ex. component = 'Core' och få ut rätt ErrorMap baserat på vad som egentligen bara är metadata
        behöver funderas på. Kommer förmodligen endast använda hyperspys plotting, får testa lite 
        
        Kommer förmodligen passa callet ner till ErrorMap där själva plottande görs från hyperspy objectet, men får se
        '''
        pass
        
class ErrorMap():
    '''
    Each ErrorMap is mostly created via the class ErrorStack 
    
    The two types of errormaps Absolute Error Map & Relative Errror Map can be found in the following variables
    .absMap contains 1 hyperspyobject which is the full navX navY|signal size
    .relMap contains one hyperspy object per element given stored in a list, basicly like the hyperspyfunction .get_lines_intesity
    
    Other metadata such as component type can be found in .type 
    
    The original data can be found in .constructed and .original
    
    '''
    def __init__(self, constructed, original, elements = ['Ag', 'Cu'] , xlines = ['Ag_La','Cu_Ka'], componentType = None ):
        
        self.type = componentType
        self.elements = elements
        self.xlines = xlines
        self.original = original
        self.constructed = setCalibration(constructed,self.original)
        self.constructed.metadata.General.title = self.type
        self.absError = SpecErrAbs2D(self.constructed, self.original)
        self.createRelElementMaps()
        self.absMap = specMapDiff(self.constructed,self.original)
        self.absMap.metadata.General.title = 'Abs. ErrorMap of: ' + self.constructed.metadata.General.title
        
        
    def __repr__(self):
        '''
        Allows for a nice representation of the object in the console
        '''
        
        return "<ErrorMap, Type: %s with Abs error %f and %s element maps>" %(self.type ,self.absError,self.elements)
    
    def __str__(self):
        '''
        Allows for a nice representation of the object in the console
        Comment:(Not sure if both are needed)
        '''
        return "<ErrorMap, Type: %s with Abs error %d and %d relative element maps>" (self.type,self.absError,len(self.elements))
    
    def createRelElementMaps(self):
        self.constructed.set_elements(self.elements)
        self.constructed.add_lines(self.xlines)
        self.original.set_elements(self.elements)
        self.original.add_lines(self.xlines)
        
        EF_maps = self.constructed.get_lines_intensity()
        refEF_maps = self.original.get_lines_intensity()
        
        self.relMaps = []
        for i in range(len(EF_maps)):
            self.relMaps.append(rel(EF_maps[i],refEF_maps[i]))
            #Constructs the metadata title of the hyperspy object 
            self.relMaps[-1].metadata.General.title = 'Rel. ErrorMap of: ' + self.constructed.metadata.General.title +' ' + EF_maps[i].metadata.Sample.xray_lines[0]

    def fetch(self,errorType = 'all', elements = ['Ag', 'Cu']):
        '''
        Default is all errortypes with Ag and Cu
        Attributes:
            errorTypes:
                'all' - self explanatory, all errortypes
                'relative' - relative errormaps divided into elements
                'absolute' - absolute errormaps only total so far
            elements - Only relevant for 'all' or relative, choose which element maps to fetch
            
            
        Fetch hyperspy objects from the errormap object based on errortype and element 
        
        
        
        Funderar på detta är nödvändigt eller om det räcker med att bara använda exempelvis ErrorStack.stack[0].relMaps osv. 
        TODO annars :)
        
        
        '''
        
        
                #if (errorType == 'all'):
        #    pass
        #pass
    
        
        #return ball
    
    def plot(self, **kwargs):
        '''
        Funderar på att ha en funktion som vid default settings plottar allt 
        annars med parametrar kan styra vilken typ av errorMap som ska plottas
        och med vilka elements (default Cu Ag typ)
        '''
        pass
    
    
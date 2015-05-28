#!/usr/bin/env python
#Classes for annotation storage information
#Dan Stribling
#Version 1.0, 01/31/2014

#from helper import *
import os.path

SPLIT_1='&&&'
SPLIT_2='!!!'
SPLIT_3='***'


# ---- Annotation Record Class ----
# Stores matches for a particular sequence from a blast result, 
#     relevant source data, and the associated annotation text.
# Stores Data:
#     ID:         String, Name of matching annotation sequence (eg. "gi|156187096|gb|EF584470.1|")
#     fileName:   String, Name of file source for match (eg. "ncbi.fasta")    
#     auxName:    String, Auxilary name (for cases of name remapping, eg. "Contig44747")
#     annotation: String, Annotation Information (eg. "Protein MGF 110-4L OS=African swine fever virus")
#     eVal:       Float, E-Score match value. (eg. 3e-52)


class AnnotationRecord():
    def __init__(self, ID, fileName=None, annotation=None, eVal=None):
        self.ID = ID
        self.fileName = fileName
        self.annotation = annotation
        self.eVal = eVal
        self.auxName = None


    # Explicitly update an existing, matching record with eVal from record provided.
    def Update(self, other):
        if not other == self:
            raise Exception('Annotation %s not equal to annotation %s.\nCannot Update' % (other.format(), self.format() ) )
        if self.eVal > other.eVal:
            self.annotation = other.annotation
            self.eVal = other.eVal

    # Return Record ID String
    def CMBID(self):
        return self.ID + ':' + self.fileName

    # Return Formatted String of E-Score Value
    def FeVal(self):
        return ('%.0e' % self.eVal).rjust(6) 

    #Return Formatted Annotation Information
    def Format(self):
        return '%s, %s, %s' % (self.ID, self.FeVal(), self.fileName)


    # --Python Magic Methods--
    # Equality Comparison, records considered equal if have the same file source and ID
    def __eq__(self, other):
        if not isinstance(other, AnnotationRecord):
            return False
        return ( (self.ID == other.ID) and (self.fileName == other.fileName) )

    # Hashing, allow to serve as dict key.
    def __hash__(self):
        return hash(self.ID + self.fileName)

    # Return as string, allows simple printing/identification.
    def __str__(self):
        return '<%s>' % self.CMBID()

    #Legacy Bridge
    update = Update
    combine_id = CMBID
    formatted_e_val = FeVal
    return_formatted = Format

#Legacy Bridge:
Annotation_Record = AnnotationRecord


# ---- Annotation Class ----
# Stores annotation records for a particular sequence in a hash table.
# The ID for each record must be unique: (Unique match name per file)
# If Non-Unique record added, will pick record with lowest e-score value.
# Stores data: 
#     name:    string, name of sequence (eg. "Contig44747")
#     records: dict, stores AnnotationRecord class instances
#         key:   string, AnnotationRecord ID (record.CMBID())
#         value: AnnotationRecord, (record)

class Annotation():
    def __init__(self, name):
        self.name = name
        self.records = {}
    
    # Add an AnnotationRecord entry to the Annotation
    # Disallows duplicate record entries (existingRecord.CMBID() == newRecord.CMBID())
    def Add(self, record):
        if record.CMBID() in self.records:
            raise Exception('Annotation ID %s already present for sequence '
                            + '%s.\nCannot Overwrite' % (record.CMBID(), self.name))
        self.records[record.CMBID()] = record

    # Update the eVal (and annotation) of an existing AnnotationRecord, if new record has smaller eVal.
    def Update(self, record):
        if record.CMBID() not in self.records:
            raise Exception('Annotation ID %s not present for sequence '
                            + '%s.\nCannot Update.' % (record.CMBID(), self.name))
        oldRecord = self.records[record.CMBID()] 
        oldRecord.Update(record)     

    # Return an AnnotationRecord by ID, if record exists in Annotation.
    def GetRecord(self, record):
        if record.CMIBD() in self.records:
            return self.records[record.CMBID()]
        else:
            return None

    # Return a list of AnnotationRecords matching subset criteria.
    # subset:    string
    #     'best', 'multibest': Returns record(s) with lowest eVal
    #         Permits multiple "best" entries with equivalent eVals
    #     'singlebest': Returns record with lowest eVal
    #         Arbitrarily chooses one entry if equivalent eVals
    # threshold: float, requires returned entries to <= threshold value.
    # sort:      boolean, return list sorted by eVal (smallest first) 
    def ReturnSubset(self, subset=None, threshold=None, sort=False):
        retList = []
        values = self.records.values()
        if subset in ['best', 'multiBest', 'singleBest']:
            if not values:
                raise Exception('No values before culling???\n%s %s' % (str(self.name), str(self.records)))
            retList =  [values[0]]
            if len(self.records) > 1:
                for record in values[1:]:
                    if (record.eVal == retList[0].eVal) and subset != 'singleBest':
                        retList.append(record)
                    elif record.eVal <= retList[0].eVal:
                        retList = [record]
        else:
            retList = values

        if threshold != None:
            if not isinstance(threshold, float):
                try:
                    threshold = float(threshold)
                except ValueError:
                    print 'Given Threshold Value (%s) Cannot Be Converted To a Float.'

            for i in reversed(range(len(retList))):
                if retList[i].eVal > threshold:
                    del(retList[i])

        if sort:
            retList = sorted(retList, key=lambda record: record.eVal)

        return retList 
 
    # Cull AnnotationRecords to those that match criteria.
    # (See ReturnSubset method)
    def Cull(self, subset=None, threshold=None):
        newRecords = {}
        for record in self.ReturnSubset(subset, threshold):
            newRecords[record.CMBID()] = record
        result = (len(self.records), len(newRecords))
        self.records = newRecords
        return result

    def Count(self, subset=None, threshold=None):
        return len(self.ReturnSubset(subset, threshold))

    def CountSources(self, subset=None, threshold=None):
        sources = {}
        #for record in self.records.values():
        for record in self.ReturnSubset(subset, threshold):
            source = "Unknown"
            if record.fileName:
                source = record.fileName

            if source in sources:
                sources[source] += 1
            else:
                sources[source] = 1

        return sources

    count_sources = CountSources
        
    # Return multi-line string with stored AnnotationRecord information
    # Indent:    integer, number of justification spaces for annotation name.
    # firstLine: boolean, put first record information on same line as annotation name.
    # subset:    See ReturnSubset method.
    # trheshold: See ReturnSubset method.
    def FormatMulti(self, indent=10, subset=None, threshold=None, firstLine=True, annotation=False):        
        out = self.name.ljust(indent)
        records = self.ReturnSubset(subset, threshold, sort=True)
        if not records:
            return ''
        if firstLine:
            out += ' : ' + records[0].Format()
            start = 1
        else:
            start = 0

        for record in records[start:]:
            out += '\n' + ''.ljust(indent) + ' : ' + record.Format()

        if not firstLine:
            out += '\n'

        return out

    # Return single-line string with stored AnnotationRecord information
    # Indent:    integer, number of justification spaces for annotation name.
    # subset:    See ReturnSubset method.
    # trheshold: See ReturnSubset method.
    def FormatSingle(self, indent=10, subset=None, threshold=None):
        printList = self.ReturnSubset(subset, threshold, sort=True)
        out = self.name.ljust(indent) + ' : '
        out += ' &!& '.join([record.Format() for record in printList])
        return out

    def FormatAnnotations(self, indent=10):
        printRecords = self.ReturnSubset(sort=True)
        if not printRecords:
            return ''
        out = ''
        for record in printRecords:
            out += self.name.ljust(indent) 
            if record.auxName:
                out += ' ' + record.auxName.ljust(indent)
            out += ' ' + record.FeVal().ljust(6) + ' ' + record.annotation + '\n'  
        return out

    def FormatFormal(self, auxName=False):
        printRecords = self.ReturnSubset(sort=True)
        out = ''
        if printRecords:
            for record in printRecords:
                out += self.name + '\t'  
                if auxName and record.auxName:
                    out += record.auxName + '\t'
                out += record.annotation + '\n'  
        return out

    def FormatList(self, subset=None, threshold=None):
        return [x.annotation for x in self.ReturnSubset(subset, threshold, sort=True)]



    # Print output of FormatMulti method (if exists)        
    def PrintFormatMulti(self, indent=10, subset=None, threshold=None, firstLine=True):
        out = self.FormatMulti(indent, subset, threshold, firstLine)
        if out:
            print out

    # Print output of FormatSingle method (if exists)
    def PrintFormatSingle(self, indent=10, subset=None, threshold=None):
        out = self.FormatSingle(indent, subset, threshold)
        if out:
            print out

    def PrintFormatAnnotations(self, indent=10):
        out = self.FormatAnnotations(indent=indent)
        if out:
            print out

    # Translate stored annotation information with specific coding.
    # Used for writing to an annotation database (annDB) file.       
    # subset:    See ReturnSubset method.
    # trheshold: See ReturnSubset method.
    def Code(self, subset=None, threshold=None):
        codeList = self.ReturnSubset(subset, threshold)
        out = self.name + SPLIT_1
        formattedList = []
        for record in codeList:
            appList = []
            for item in [record.ID, str(record.eVal), record.fileName, record.annotation]:
                if item not in [None, '']:
                    appList.append(item)
 
            formattedList.append(SPLIT_2.join(appList))

        out += SPLIT_3.join(formattedList)
        return out

    # Translate coded annotation information into class instance.
    # Used for reading from an annotation database (annDB) file.       
    def Decode(self, line):
        lineSplit=line.split(SPLIT_1)
        annotation = Annotation(name=lineSplit[0])
       
        rawRecords = lineSplit[1].split(SPLIT_3)

        for rawRecord in rawRecords:
            rawSplit = rawRecord.split(SPLIT_2)
            #eprint rawSplit
            record=AnnotationRecord(ID=rawSplit[0])
            if len(rawSplit) > 1:
                record.eVal = float(rawSplit[1])
            if len(rawSplit) > 2:
                record.fileName = rawSplit[2]
            if len(rawSplit) > 3:
                record.annotation = rawSplit[3]
            annotation.Add(record) 
        return annotation      


    # --Python Magic Methods--
    # Enable iteration over class instance  
    def __iter__(self):
        return self.records.values().__iter__()    

    # Equality comparison, based on whether Annotation instances have same name.
    def __eq__(self, other):
        if not isistance(other, Annotation):
            return False
        return self.name == other.name

    #def __hash__(self):
    #    return hash(self.name)

    # Return String of Annotation instance
    def __str__(self):
        return '<ANNOTATION:%s>' % self.name

    # Length datum, return nubmer of records in Annotation
    def __len__(self):
        return len(self.records)

    # Boolean comparison, returns True if instance contains records not empty.
    def __nonzero__(self):
        return bool(self.records)

    #Legacy Bridge
    add = Add
    update = Update
    get_record = GetRecord
    return_subset = ReturnSubset
    cull = Cull
    count = Count
    count_sources = CountSources
    format_multi = FormatMulti
    format_single = FormatSingle
    format_annotations = FormatAnnotations
    format_formal = FormatFormal
    format_list = FormatList
    print_format_multi = PrintFormatMulti
    print_format_single = PrintFormatSingle
    print_format_annotations = PrintFormatAnnotations
    code = Code
    decode = Decode





# ---- Annotation Database ----
# Stores Annotation instances in a hash table
# Requires unique annotation names for correct functioning.
# Stores data: 
#     fileName:    string, name of last database file read/write
#     annotations: dict, stores Annotation class instances
#         key:   string, Annotation name (annotation.name)
#         value: Annotation, (annotation)

class AnnotationDB():
    def __init__(self, fileName = None):
        self.annotations = {}
        self.fileName = None
        if fileName:
            self.Read(fileName)

    # Read information from coded annotation database file. (*.annDB)
    def Read(self, fileName):
        if not os.path.isfile(fileName):
            raise Exception('Problem Database Read File %s Not Found' % fileName)
        inFile = open(fileName, 'r')
        for line in inFile:
            self.AddAnnotation(Annotation(None).Decode(line.rstrip())) 
        inFile.close()

    # Write information to coded annotation database file. (*.annDB) 
    def Write(self, fileName=None, subset=None):
        if not fileName:
            if not self.fileName:
                raise Exception('No output filename provided for writing.')
            fileName = self.fileName

        outFile = open(fileName, 'w')
        for annotation in self._SortedAnnotations():
            #print annotation.Code()
            outFile.write(annotation.Code(subset=subset) + '\n')
        outFile.close()

    # Add Annotation class instance to database.
    # Disallowed addition of duplicates.        
    def AddAnnotation(self, annotation):
        if not isinstance(annotation, Annotation):
            raise Exception('Only annotations can be added to AnnotationDB via Add Method')
        if annotation.name in self.annotations:
            raise Exception('Annotation: %s already exists in database.'
                            + 'Unique annotation names required.' % annotation.name)
        self.annotations[annotation.name]=annotation

    # Access Annotation class instance from database by name and return.
    def GetAnnotation(self, name):
        if name not in self:
            raise Exception('Annotation %s not found in database.' % name)
        return self.annotations[name]

    # Adds AnnotationRecord to appropriate Annotation instance in database. 
    # Creates Annotation instance in Database if necessary.
    # Allowes duplicate record entries within annotations by using Annotation.Update
    #     method if AnnotationRecord already exists in Annotation instance.
    def AddRecord(self, name, record):
        if not isinstance(name, str):
            raise Exception('name: %s not a string' % str(name))
        if not isinstance(record, AnnotationRecord):
            raise Exception('record: %s is not type AnnotationRecord' % str(record))      

        if name not in self.annotations:
            self.annotations[name] = Annotation(name=name)
            self.annotations[name].Add(record)
        else:
            if record in self.annotations[name]:
                self.annotations[name].Update(record)
            else:
                self.annotations[name].Add(record)

    # Not Implemented
    #def _GetRecord(self, name, recordID):
    #    if not isinstance(name, str):
    #        raise Exception('name: %s not a string' % str(name))
    #    if not isinstance(recordID, str):
    #        raise Exception('recordID: %s is not a string' % str(record))      
    #    if name not in self.annotations:
    #        raise Exception('Annotation: %s not found' % str(name))


    # Cull all annotations in AnnotationDB to specified criteria.
    # Returns tuple: (
    #     initial number of Annotation instances,
    #     initial total number of AnnotationRecord instances
    #     final number of Annotation instances,
    #     final total number of AnnotationRecord instances         
    # subset:    See ReturnSubset method of Annotation class.
    # trheshold: See ReturnSubset method of Annotation class.
    def Cull(self, subset=None, threshold=None):
        initialAnnotations = len(self.annotations)
        initialRecords = 0
        finalRecords = 0
        if subset or threshold != None:
            for key in self.annotations.keys():
                (ini, fin) = self.annotations[key].Cull(subset, threshold)
                initialRecords += ini
                finalRecords += fin
                if not fin:
                    del(self.annotations[key])

        finalAnnotations = len(self.annotations)

        return (initialAnnotations, initialRecords, finalAnnotations, finalRecords) 

    def Count(self, subset=None, threshold=None):
        records = 0
        if subset or (threshold != None):
            for annotation in self.annotations.values():
                records += annotation.Count(subset, threshold)
        else:
            for annotation in self.annotations.values():
                records += len(annotation)

        return (len(self.annotations), records)             


    def CountSources(self, subset=None, threshold=None, one_count=True):
        #sources = {'total':0}
        sources = {}
        for annotation in self.annotations.values():
            annotation_sources = annotation.CountSources(subset=subset, threshold=threshold)
            for source in annotation_sources:
                if source in sources:
                    if one_count:
                        sources[source] += 1
                    else:
                        sources[source] += annotation_sources[source]
                else:
                    if one_count:
                        sources[source] = 1
                    else:
                        sources[source] = annotation_sources[source]
                #sources['total'] += annotation_sources[source]

        return sources



    # Return list of stored Annotation instances, sorted by Annotation name
    # Use natsort library for "natural" sorting if available, else use "standard" sorting. 
    def _SortedAnnotations(self):
        try:
            import natsort
            return natsort.natsorted(self.annotations.values(), key=lambda annotation: annotation.name) 
        except ImportError:
            return sorted(self.annotations.values(), key=lambda annotation: annotation.name)
 
    # ??? Found 03-23-2015
    #def _ReturnSubset(self, subset):
    #    aa


    # Uses AnnotationMap class instance to add annotation strings to Annotation instances in databse.
    # Expects to find an annotation for each instance, or will raise error.
    def MapAnnotations(self, annotationMap):
        if not isinstance(annotationMap, AnnotationMap):
            raise Exception('Annotation map provided: %s' % str(annotationMap)
                            + ' is not an instance of AnnotationMap Class.') 
        additions = 0
        for annotation in self.annotations.values():
            for record in annotation:
                record.annotation = annotationMap[record.ID]
                additions += 1
        return additions
   
    # Uses NameMap class instance to remap annotations onto new names.
    # If no new name provided, will use old name.
    # If two annotations are mapped onto the same new name, they will be merged.
    def MapNames(self, nameMap, debug=False):
        if not isinstance(nameMap, NameMap):
            raise Exception('Name map provided: %s' % str(nameMap)
                            + ' is not an instance of NameMap Class.') 
        newAnnotations = {}
        startingAnnotations = len(self.annotations)
        for annotation in self.annotations.values():
            oldName = annotation.name

            for record in annotation:
                record.auxName = oldName

            if oldName in nameMap:
                newName = nameMap[oldName]

                #if newName == 'Gene040810':
                    #print newName
                    #print oldName
                    #print annotation
                    #print annotation.records
                    #raw_input()


                if newName in newAnnotations:
                    for record in annotation:
                        if record in newAnnotations[newName]:
                            newAnnotations[newName].Update(record)
                        else:
                            newAnnotations[newName].Add(record)
                    if debug:
                        print 'Adding %s to existing annotation: %s' % (oldName, newName)
             
                else:
                    annotation.name = newName
                    newAnnotations[newName] = annotation
                    if debug:
                        print 'Mapping %s to %s' % (oldName, newName)
            else:
                newAnnotations[oldName] = annotation

        self.annotations = newAnnotations 
       
        endingAnnotations = len(self.annotations)
        return (startingAnnotations, endingAnnotations)    

    # Combine two AnnotationDB instances.
    def Combine(self, otherDB):
        if not isinstance(otherDB, AnnotationDB):
            raise Exception('AnnotationDB provided: %s' % str(otherDB)
                            + ' is not an instance of AnnotationDB Class.') 
        for annotation in otherDB:
            if annotation.name not in self.annotations:
                self.annotations[annotation.name] = annotation
            else:
                for record in annotation:
                    self.AddRecord(annotation.name, record)
        
 

    # Return string with sorted AnnotationRecord information stored in database
    # Indent:     integer, number of justification spaces for annotation name.
    # subset:     See ReturnSubset method of Annotation class.
    # trheshold:  See ReturnSubset method of Annotation class.
    # firstLine:  boolean, put first record information on same line as annotation name.
    # formatType: string, line arrangement for Annotations 
    # spacing:    boolean, add one line of space between Annotations
    def Format(self, indent=10, subset=None, threshold=None, firstLine=True, formatType='multi', spacing=False):
        retStr=''

        if formatType == 'multi':
            for annotation in self._SortedAnnotations():
                annStr = annotation.FormatMulti(indent, subset, threshold, firstLine)
                if annStr:
                    retStr += annStr + '\n'
                    if spacing:
                        retStr += '\n'

        elif formatType == 'single':
            for annotation in self._SortedAnnotations():
                annStr = annotation.FormatSingle(indent, subset, threshold)
                if annStr:
                    retStr += annStr + '\n'
                    if spacing:
                        retStr += '\n'
        
        elif formatType == 'annotation':
            for annotation in self._SortedAnnotations():
                annStr = annotation.FormatAnnotations(indent)
                if annStr:
                    retStr += annStr
                    if spacing:
                        retStr += '\n'

        elif formatType == 'formal':
            for annotation in self._SortedAnnotations():
                annStr = annotation.FormatFormal()
                if annStr:
                    retStr += annStr
                    if spacing:
                        retStr += '\n'



        else:
            raise Exception('Unrecognized AnnotationDB Output formatting type: %s' % str(formatType))

        return retStr

    # Print output of Format method, if it exists.
    def PrintFormat(self, indent=10, subset=None, threshold=None, firstLine=True, formatType='multi', spacing=False):
        out = self.Format(indent, subset, threshold, firstLine, formatType, spacing)
        if out:
            print out
        else:
            print 'No Annotations Match Criteria'


    # --Python Magic Methods--
    # Item in container:
    def __contains__(self, name):
        return name in self.annotations        
    
    # Number of Annotation class instances in database.
    def __len__(self):
        return len(self.annotations) 

    # Iterate over Annotation class instances in database.
    def __iter__(self):
        return iter(self.annotations.values())

    #Legacy Bridge
    read = Read
    write = Write
    add_annotation = AddAnnotation
    get_annotation = GetAnnotation
    add_record = AddRecord
    cull = Cull
    count = Count
    count_sources = CountSources
    map_annotations = MapAnnotations
    map_names = MapNames
    combine = Combine
    return_formatted = Format
    print_format = PrintFormat

Annotation_Database = AnnotationDB



# ---- AnnotationMap Class ----
# Stores sequence headers from a sequence file.
#     ID:         String, Name of matching annotation sequence (eg. "gi|156187096|gb|EF584470.1|")
#     fileName:   String, Name of file source for match (eg. "ncbi.fasta")    
#     annotation: String, Annotation Information (eg. "Protein MGF 110-4L OS=African swine fever virus")
class AnnotationMap():
    def __init__(self, fileName=None):  
        self.annotations = {}
        self.fileName = None
        if fileName:
            self.Read(fileName)

    # Read information from sequence file. (*.annDB)
    def Read(self, fileName):
        if not os.path.isfile(fileName):
            raise Exception('Problem: Annotation Map Read File %s Not Found' % fileName)
        if fileName.lower().endswith(('.fa', '.fas', '.fasta', '.FA', '.FASTA')):
            self._ReadFASTA(fileName)
        else:
            raise Exception('Type for file: %s not supported.' % fileName
                             + '\n(Be sure file ends with .fa, .fasta, etc.')
        self.fileName = fileName


    def _ReadFASTA(self, fileName):
        inFile = open(fileName, 'r')
        for line in inFile:
            if line.startswith('>'):
                splitLine = line.lstrip('>').split()
                if len(splitLine) < 2:
                    raise Exception('Line has insufficent information:\n' + line)
                self.annotations[splitLine[0]] = ' '.join(splitLine[1:])                            
        inFile.close()
    
    # --Python Magic Methods--
    # Item in container:
    def __contains__(self, name):
        return name in self.annotations        
    
    # Number of Annotation class instances in database.
    def __len__(self):
        return len(self.annotations) 

    # Iterate over Annotation class instances in database.
    def __iter__(self):
        return iter(self.annotations.values())

    # Return whether AnnotationMap is empty.
    def __nonzero__(self):
        return bool(self.annotations)

    # Retrieve annotation by name.
    def __getitem__(self, name):
        return self.annotations[name]

    #Legacy Bridge
    read = Read

Annotation_Map = AnnotationMap





# ---- NameMap Class ----
# Used for renaming sequence headers from a sequence file.
# Stores startName:endName information in dictionary.
# Takes a custom parser function as an initizliation arguement.
# Parse function must take a string and return a (startName, endName) tuple. 
#     startName: String, Starting map name. (eg. "Contig44747")
#     startName: String, Ending map name. (eg. "Gene00072")
#     fileName:  String, Name of file source for match (eg. "ncbi.fasta")    

# Default Parser function for Name Map Class
def DefaultMapParse(inputStr):
    inputStr = inputStr.strip()
    splitString = inputStr.split('|')
    return (splitString[0].lstrip('>').strip(), splitString[1].strip())
default_map_parse = DefaultMapParse

class NameMap():
    def __init__(self, fileName=None, parser=DefaultMapParse):  
        self.nameMap = {}
        self.fileName = None
        self.parser = parser
        if fileName:
            self.Read(fileName)

    # Read information from file.)
    def Read(self, fileName):
        if not os.path.isfile(fileName):
            raise Exception('Problem: Name Map Read File %s Not Found' % fileName)
        if fileName.lower().endswith(('.fa', '.fasta', '.FA', '.FASTA')):
            self._ReadFASTA(fileName)
        else:
            self._ReadGeneric(fileName)
        self.fileName = fileName
        return len(self.nameMap)

    def _ReadFASTA(self, fileName):
        inFile = open(fileName, 'r')
        for line in inFile:
            if line.startswith('>'):
                (startName, endName) = self.parser(line)
                self.nameMap[startName] = endName                           
        inFile.close()
    
    def _ReadGeneric(self, fileName):
        inFile = open(fileName, 'r')
        for line in inFile:
            (startName, endName) = self.parser(line)
            self.nameMap[startName] = endName                           
        inFile.close()
       

    # --Python Magic Methods--
    # Item in container:
    def __contains__(self, name):
        return name in self.nameMap
    
    # Number of names mapped in instance.
    def __len__(self):
        return len(self.nameMap) 

    # Iterate over name maps in instance.
    def __iter__(self):
        return iter(self.nameMap)

    # Return whether instance is empty.
    def __nonzero__(self):
        return bool(self.nameMap)

    # Retrieve instance by name.
    def __getitem__(self, name):
        return self.nameMap[name]

    #Legacy Bridge
    read = Read

Name_Map = NameMap

def FindSequenceAnnotations(targetSequences, blastResults, annotationDBName):

    if not os.path.isfile(targetSequences):
        print '\nTarget Sequence File %s Not Found.\n' % fileName
        sys.exit(1)

    if not os.path.isfile(annotationDBName):
        print '\nAnnotation Database File %s Not Found.\n' % annotationDBName
        sys.exit(1)

    if not os.path.isfile(blastResults):
        print '\nBlast Results File: %s Not Found.\n' % blastResults
        sys.exit(1)
   
    print 'Finidng Annotations for File: %s' % targetSequences 
    print 'Using Annotation Database: %s' % annotationDBName
    print 'Using BLAST output file: %s' % blastResults

    print '\nReading Annotations from Database:' 

    annotationDB = {}
    annotationFile = open(annotationDBName, 'r')
    for line in annotationFile:
        splitLine = line.split('-!-')
        #print line
        #print splitLine
        annotationDB[splitLine[0]] = splitLine[1]
    annotationFile.close()

    matchDB = {}
    resultsFile = open(blastResults, 'r')
    for line in resultsFile:
        splitLine = line.split()
        gene = splitLine[0]
        sequence = splitLine[1]
        eScore = float(splitLine[10])
        if gene not in annotationDB:
            print '\nPROBLEM!!! Gene %s Not Found in Annotation Database: %s' % (gene, annotationDB)
            print splitLine
            sys.exit(1)

        if eScore <= ESCORE_CUTOFF:
            if sequence in matchDB:
               print 'Gene: %s matches previously annotated sequence %s, Skipping.' % (gene, sequence)
             
            matchDB[sequence] = annotationDB[gene]
 
    resultsFile.close()

find_sequence_annotations = FindSequenceAnnotations

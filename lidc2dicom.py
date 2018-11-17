from pathlib import Path
import lidc_conversion_utils.helpers as lidc_helpers
import os, itk, tempfile, json, pydicom, tempfile, shutil
from subprocess import call
import pylidc as pl
import numpy as np
import glob
import logging
from decimal import *

class LIDC2DICOMConverter:

  def __init__(self):
    self.logger = logging.getLogger("lidc2dicom")

    self.rootDir = "/Users/fedorov/Documents/TCIA/LIDC-IDRI"
    self.segTemplate = "seg_conversion_template.json"
    self.srTemplate = "sr_conversion_template.json"
    self.colorsFile = "GenericColors.txt"
    self.tempDir="/Users/fedorov/Temp/LIDC_conversion1"
    self.dcmqiRoot = "/Users/fedorov/local/dcmqi-1.1.0-mac-20181017-ef688c4/bin"

    # read GenericColors
    self.colors = []
    with open(self.colorsFile,'r') as f:
      for l in f:
        if l.startswith('#'):
          continue
        self.colors.append([int(c) for c in l.split(' ')[2:5]])

  def cleanUpTempDir(self, dir):
    shutil.rmtree(os.path.join(dir,"dicom2nrrd"))
    for p in Path(dir).glob("*.nrrd"):
      p.unlink()

  def convertSingleAnnotation(self, nCount, aCount, a, ctDCM, noduleUID, volume, seriesDir):
    with open(self.segTemplate,'r') as f:
      segJSON = json.load(f)

    # update as necessary!
    segName = "Nodule "+str(nCount+1) +" - Annotation " + a._nodule_id
    segJSON["segmentAttributes"][0][0]["SegmentDescription"] = segName
    segJSON["SeriesDescription"] = "Segmentation of "+segName

    if ctDCM.SeriesNumber != '':
      segJSON["SeriesNumber"] = str(int(ctDCM.SeriesNumber)+self.instanceCount+1)
    else:
      segJSON["SeriesNumber"] = str(self.instanceCount+1)

    for ci in range(3):
        segJSON["segmentAttributes"][0][0]["recommendedDisplayRGBValue"][ci] = self.colors[aCount+1][ci]

    segJSON["segmentAttributes"][0][0]["TrackingIdentifier"] = segName
    segJSON["segmentAttributes"][0][0]["TrackingUniqueIdentifier"] = noduleUID

    jsonSegFile = os.path.join(self.tempSubjectDir,segName+'.json')
    with open(jsonSegFile, "w") as f:
      json.dump(segJSON, f, indent=2)

    maskArray = a.boolean_mask(10000).astype(np.int16)

    maskArray = np.swapaxes(maskArray,0,2).copy()
    maskArray = np.rollaxis(maskArray,2,1).copy()

    maskVolume = itk.GetImageFromArray(maskArray)
    maskVolume.SetSpacing(volume.GetSpacing())
    maskVolume.SetOrigin(volume.GetOrigin())
    writerType = itk.ImageFileWriter[itk.Image[itk.SS, 3]]
    writer = writerType.New()
    nrrdSegFile = os.path.join(self.tempSubjectDir,segName+'.nrrd')
    writer.SetFileName(nrrdSegFile)
    writer.SetInput(maskVolume)
    writer.Update()

    dcmSegFile = os.path.join(self.tempSubjectDir,segName+'.dcm')

    converterCmd = [os.path.join(self.dcmqiRoot,'itkimage2segimage'), "--inputImageList", nrrdSegFile, "--inputDICOMDirectory", seriesDir, "--inputMetadata", jsonSegFile, "--outputDICOM", dcmSegFile]
    self.logger.info("Converting to DICOM SEG with "+str(converterCmd))
    call(converterCmd)

    segUID = None
    ctSeriesUID = None
    try:
      segDcm = pydicom.read_file(dcmSegFile)
      segUID = segDcm.SOPInstanceUID
      ctSeriesUID = segDcm.ReferencedSeriesSequence[0].SeriesInstanceUID
    except:
      self.logger.error("Failed to read Segmentation file")
      return

    self.instanceCount = self.instanceCount+1

    # convert qualitative evaluations
    # TODO: DAC suggests it is better mix and match standard codes than use non-standard for all
    lidcCodingScheme = "99LIDCQIICR"
    conceptsDictionary = {'subtlety':{"CodeValue":'C45992',"CodeMeaning":'Subtlety score',"CodingSchemeDesignator":"NCIt"},
                          'internalStructure':{"CodeValue":'200',"CodeMeaning":'Internal structure',"CodingSchemeDesignator":lidcCodingScheme},
                          'calcification':{"CodeValue":'C3672',"CodeMeaning":"Calcification","CodingSchemeDesignator":"NCIt"},
                          'sphericity':{"CodeValue":'400',"CodeMeaning":"Sphericity","CodingSchemeDesignator":lidcCodingScheme},
                          'margin':{"CodeValue":'C25563',"CodeMeaning":'Margin',"CodingSchemeDesignator":"NCIt"},
                          'lobulation':{"CodeValue":'600',"CodeMeaning":'Lobulation',"CodingSchemeDesignator":lidcCodingScheme},
                          'spiculation':{"CodeValue":'700',"CodeMeaning":'Spiculation',"CodingSchemeDesignator":lidcCodingScheme},
                          'texture':{"CodeValue":'C41144',"CodeMeaning":'Texture',"CodingSchemeDesignator":"NCIt"},
                          'malignancy':{"CodeValue":'900',"CodeMeaning":'Malignancy',"CodingSchemeDesignator":lidcCodingScheme}}

    valuesDictionary = {}

    valuesDictionary['generic'] = {}
    valuesDictionary['generic']['1'] = {"CodeValue":'001',"CodeMeaning":"1 out of 5","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['generic']['2'] = {"CodeValue":'002',"CodeMeaning":"2 out of 5","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['generic']['3'] = {"CodeValue":'003',"CodeMeaning":"3 out of 5","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['generic']['4'] = {"CodeValue":'004',"CodeMeaning":"4 out of 5","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['generic']['5'] = {"CodeValue":'005',"CodeMeaning":"5 out of 5","CodingSchemeDesignator":lidcCodingScheme}


    valuesDictionary['subtlety'] = {}
    valuesDictionary['subtlety']['1'] = {"CodeValue":'101',"CodeMeaning":"1 out of 5 (Extremely subtle)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['subtlety']['2'] = {"CodeValue":'102',"CodeMeaning":"2 out of 5 (Moderately subtle)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['subtlety']['3'] = {"CodeValue":'103',"CodeMeaning":"3 out of 5 (Fairly subtle)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['subtlety']['4'] = {"CodeValue":'104',"CodeMeaning":"4 out of 5 (Moderately obvious)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['subtlety']['5'] = {"CodeValue":'105',"CodeMeaning":"5 out of 5 (Obvious)","CodingSchemeDesignator":lidcCodingScheme}

    valuesDictionary['internalStructure'] = {}
    valuesDictionary['internalStructure']['1'] = {"CodeValue":'C12471',"CodeMeaning":"Soft tissue","CodingSchemeDesignator":"NCIt"}
    valuesDictionary['internalStructure']['2'] = {"CodeValue":'C25278',"CodeMeaning":"Fluid","CodingSchemeDesignator":"NCIt"}
    valuesDictionary['internalStructure']['3'] = {"CodeValue":'C12472',"CodeMeaning":"Adipose tissue","CodingSchemeDesignator":"NCIt"}
    valuesDictionary['internalStructure']['4'] = {"CodeValue":'C73434',"CodeMeaning":"Air","CodingSchemeDesignator":"NCIt"}

    # TODO: revisit RadLex and DCM for those codes
    valuesDictionary['calcification'] = {}
    valuesDictionary['calcification']['1'] = {"CodeValue":'RID35453',"CodeMeaning":"Popcorn calcification sign","CodingSchemeDesignator":"RadLex"}
    valuesDictionary['calcification']['2'] = {"CodeValue":'302',"CodeMeaning":"Laminated appearance","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['calcification']['3'] = {"CodeValue":'303',"CodeMeaning":"Solid appearance","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['calcification']['4'] = {"CodeValue":'304',"CodeMeaning":"Non-central appearance","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['calcification']['5'] = {"CodeValue":'305',"CodeMeaning":"Central calcification","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['calcification']['6'] = {"CodeValue":'RID28473',"CodeMeaning":"Absent","CodingSchemeDesignator":"RadLex"}

    valuesDictionary['sphericity'] = {}
    valuesDictionary['sphericity']['1'] = {"CodeValue":'RID5811',"CodeMeaning":"linear","CodingSchemeDesignator":"RadLex"}
    valuesDictionary['sphericity']['2'] = valuesDictionary['generic']['2']
    valuesDictionary['sphericity']['3'] = {"CodeValue":'RID5800',"CodeMeaning":"ovoid","CodingSchemeDesignator":"RadLex"}
    valuesDictionary['sphericity']['4'] = valuesDictionary['generic']['4']
    valuesDictionary['sphericity']['5'] = {"CodeValue":'RID5799',"CodeMeaning":"round","CodingSchemeDesignator":"RadLex"}

    valuesDictionary['margin'] = {}
    valuesDictionary['margin']['1'] = {"CodeValue":'RID5709',"CodeMeaning":"Indistinct margin","CodingSchemeDesignator":"RadLex"}
    valuesDictionary['margin']['2'] = valuesDictionary['generic']['2']
    valuesDictionary['margin']['3'] = valuesDictionary['generic']['3']
    valuesDictionary['margin']['4'] = valuesDictionary['generic']['4']
    valuesDictionary['margin']['5'] = {"CodeValue":'RID5707',"CodeMeaning":"Circumscribed margin","CodingSchemeDesignator":"RadLex"}

    valuesDictionary['lobulation'] = {}
    valuesDictionary['lobulation']['1'] = {"CodeValue":'601',"CodeMeaning":"1 out of 5 (No lobulation)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['lobulation']['2'] = valuesDictionary['generic']['2']
    valuesDictionary['lobulation']['3'] = valuesDictionary['generic']['3']
    valuesDictionary['lobulation']['4'] = valuesDictionary['generic']['4']
    valuesDictionary['lobulation']['5'] = {"CodeValue":'605',"CodeMeaning":"5 out of 5 (Marked lobulation)","CodingSchemeDesignator":lidcCodingScheme}

    valuesDictionary['spiculation'] = {}
    valuesDictionary['spiculation']['1'] = {"CodeValue":'701',"CodeMeaning":"1 out of 5 (No spiculation)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['spiculation']['2'] = valuesDictionary['generic']['2']
    valuesDictionary['spiculation']['3'] = valuesDictionary['generic']['3']
    valuesDictionary['spiculation']['4'] = valuesDictionary['generic']['4']
    valuesDictionary['spiculation']['5'] = {"CodeValue":'705',"CodeMeaning":"5 out of 5 (Marked spiculation)","CodingSchemeDesignator":lidcCodingScheme}

    valuesDictionary['texture'] = {}
    valuesDictionary['texture']['1'] = {"CodeValue":'RID50153',"CodeMeaning":"non-solid pulmonary nodule","CodingSchemeDesignator":"RadLex"}
    valuesDictionary['texture']['2'] = valuesDictionary['generic']['2']
    valuesDictionary['texture']['3'] = {"CodeValue":'RID50152',"CodeMeaning":"part-solid pulmonary nodule","CodingSchemeDesignator":"RadLex"}
    valuesDictionary['texture']['4'] = valuesDictionary['generic']['4']
    valuesDictionary['texture']['5'] = {"CodeValue":'RID50151',"CodeMeaning":"solid pulmonary nodule","CodingSchemeDesignator":"RadLex"}

    valuesDictionary['malignancy'] = {}
    valuesDictionary['malignancy']['1'] = {"CodeValue":'901',"CodeMeaning":"1 out of 5 (Highly Unlikely for Cancer)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['malignancy']['2'] = {"CodeValue":'902',"CodeMeaning":"2 out of 5 (Moderately Unlikely for Cancer)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['malignancy']['3'] = {"CodeValue":'903',"CodeMeaning":"3 out of 5 (Indeterminate Likelihood)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['malignancy']['4'] = {"CodeValue":'904',"CodeMeaning":"4 out of 5 (Moderately Suspicious for Cancer)","CodingSchemeDesignator":lidcCodingScheme}
    valuesDictionary['malignancy']['5'] = {"CodeValue":'905',"CodeMeaning":"5 out of 5 (Highly Suspicious for Cancer)","CodingSchemeDesignator":lidcCodingScheme}

    with open(self.srTemplate,'r') as f:
      srJSON = json.load(f)

    srName = segName+" evaluations"
    srJSON["SeriesDescription"] = srName

    srJSON["compositeContext"] = [dcmSegFile.split('/')[-1]]
    srJSON["imageLibrary"] = os.listdir(seriesDir)

    qualitativeEvaluations = []
    measurementItems = []

    volumeItem = {}
    volumeItem["value"] = '%E' % Decimal(a.volume)
    volumeItem["quantity"] = {"CodeValue": "G-D705","CodingSchemeDesignator": "SRT","CodeMeaning": "Volume"}
    volumeItem["units"] = {"CodeValue": "mm3","CodingSchemeDesignator": "UCUM","CodeMeaning": "cubic millimeter"}
    volumeItem["measurementModifier"] = {"CodeValue": "122503","CodingSchemeDesignator": "DCM","CodeMeaning": "Integration of sum of closed areas on contiguous slices"}
    volumeItem["measurementAlgorithmIdentification"] = {"AlgorithmName": "pylidc","AlgorithmVersion": "0.2.0"}

    # CID 7470
    diameterItem = {}
    diameterItem["value"] = '%E' % Decimal(a.diameter)
    diameterItem["quantity"] = {"CodeValue": "M-02550","CodingSchemeDesignator": "SRT","CodeMeaning": "Diameter"}
    diameterItem["units"] = {"CodeValue": "mm","CodingSchemeDesignator": "UCUM","CodeMeaning": "millimeter"}
    diameterItem["measurementAlgorithmIdentification"] = {"AlgorithmName": "pylidc","AlgorithmVersion": "0.2.0"}

    #
    surfaceItem = {}
    surfaceItem["value"] = '%E' % Decimal(a.surface_area)
    surfaceItem["quantity"] = {"CodeValue": "C0JK","CodingSchemeDesignator": "IBSI","CodeMeaning": "Surface area of mesh"}
    surfaceItem["units"] = {"CodeValue": "mm2","CodingSchemeDesignator": "UCUM","CodeMeaning": "square millimeter"}
    surfaceItem["measurementAlgorithmIdentification"] = {"AlgorithmName": "pylidc","AlgorithmVersion": "0.2.0"}

    measurementItems.append(volumeItem)
    measurementItems.append(diameterItem)
    measurementItems.append(surfaceItem)

    #      subtletyItem = {}
    #      subtletyItem["conceptCode"] = conceptsDictionary['subtlety']
    #      subtletyItem["conceptValue"] = valuesDictionary['subtlety'][str(a.subtlety)]
    #      qualitativeEvaluations.append(subtletyItem)

    for attribute in conceptsDictionary.keys():
      # print(attribute+': '+str(getattr(a, attribute)))

      try:
        qItem = {}
        qItem["conceptCode"] = conceptsDictionary[attribute]
        qItem["conceptValue"] = valuesDictionary[attribute][str(getattr(a, attribute))]
        qualitativeEvaluations.append(qItem)
      except KeyError:
        self.logger.info("Skipping invalid attribute: "+attribute+': '+str(getattr(a, attribute)))
        continue

    srJSON["Measurements"][0]["measurementItems"] = measurementItems
    srJSON["Measurements"][0]["qualitativeEvaluations"] = qualitativeEvaluations
    srJSON["Measurements"][0]["segmentationSOPInstanceUID"] = segUID
    srJSON["Measurements"][0]["SourceSeriesForImageSegmentation"] = ctSeriesUID

    srJSON["Measurements"][0]["TrackingIdentifier"] = segName
    srJSON["Measurements"][0]["TrackingUniqueIdentifier"] = noduleUID

    srName = "Nodule "+str(nCount+1) +" - Annotation " + a._nodule_id + " measurements"
    jsonSRFile = os.path.join(self.tempSubjectDir,srName+'.json')
    with open(jsonSRFile, "w") as f:
      json.dump(srJSON, f, indent=2)

    dcmSRFile = os.path.join(self.tempSubjectDir,srName+'.dcm')
    converterCmd = [os.path.join(self.dcmqiRoot,'tid1500writer'), "--inputMetadata", jsonSRFile, "--inputImageLibraryDirectory", seriesDir, "--inputCompositeContextDirectory", self.tempSubjectDir, "--outputDICOM", dcmSRFile]
    self.logger.info("Converting with "+str(converterCmd))
    call(converterCmd)

    if not os.path.exists(dcmSRFile):
      self.logger.error("Failed to access output SR file for "+s)

  def convertForSubject(self, subjectID):
    s = 'LIDC-IDRI-%04i' % subjectID
    self.logger.info("Processing subject %s" % (s))
    scans = pl.query(pl.Scan).filter(pl.Scan.patient_id == s)
    self.logger.info(" Found %d scans" % (scans.count()))

    for scan in scans:
      studyUID = scan.study_instance_uid
      seriesUID = scan.series_instance_uid
      seriesDir = os.path.join(self.rootDir,s,studyUID,seriesUID)
      if not os.path.exists(seriesDir):
        self.logger.error("Files not found for subject "+s)
        return

      dcmFiles = glob.glob(os.path.join(seriesDir,"*.dcm"))
      if not len(dcmFiles):
        logger.error("No DICOM files found for subject "+s)
        return

      firstFile = os.path.join(seriesDir,dcmFiles[0])

      try:
        ctDCM = pydicom.read_file(firstFile)
      except:
        logger.error("Failed to read input file "+firstFile)
        return

      ok = lidc_helpers.checkSeriesGeometry(seriesDir)
      if not ok:
        self.logger.warning("Geometry inconsistent for subject %s" % (s))

      self.tempSubjectDir = os.path.join(self.tempDir,s)
      reconTempDir = os.path.join(self.tempSubjectDir,"dicom2nrrd")
      try:
        os.makedirs(reconTempDir)
      except:
        pass

      scanNRRDFile = os.path.join(self.tempSubjectDir,s+'_CT.nrrd')
      if not os.path.exists(scanNRRDFile):
        # convert
        # tempDir = tempfile.mkdtemp()
        plastimatchCmd = ['/Users/fedorov/build/plastimatch/plastimatch', 'convert','--input',seriesDir,'--output-img',scanNRRDFile]
        self.logger.info("Running plastimatch with "+str(plastimatchCmd))
        call(plastimatchCmd)
        self.logger.info('plastimatch completed')
        self.logger.info("Conversion of CT volume OK - result in "+scanNRRDFile)
      else:
        self.logger.info(scanNRRDFile+" exists. Not rerunning volume reconstruction.")

      reader = itk.ImageFileReader[itk.Image[itk.SS, 3]].New()
      reader.SetFileName(scanNRRDFile)
      reader.Update()
      volume = reader.GetOutput()

      #logger.info(volume.GetLargestPossibleRegion().GetSize())

      # now iterate over all nodules available for this subject
      anns = scan.annotations
      self.logger.info("Have %d annotations for subject %s" % (len(anns), s))

      self.instanceCount = 0

      clusteredAnnotationIDs = []

      for nCount,nodule in enumerate(scan.cluster_annotations()):

        noduleUID = pydicom.uid.generate_uid(prefix=None) # by default, pydicom uses 2.25 root

        for aCount,a in enumerate(nodule):

          clusteredAnnotationIDs.append(a.id)
          self.convertSingleAnnotation(nCount, aCount, a, ctDCM, noduleUID, volume, seriesDir)


      if len(clusteredAnnotationIDs) != len(anns):
        self.logger.warning("%d annotations unaccounted for!" % (len(anns) - len(clusteredAnnotationIDs)))

      for ua in anns:
        if ua.id not in clusteredAnnotationIDs:
          aCount = aCount+1
          nCount = nCount+1
          noduleUID = pydicom.uid.generate_uid(prefix=None)
          self.convertSingleAnnotation(nCount, aCount, ua, ctDCM, noduleUID, volume, seriesDir)

      self.cleanUpTempDir(self.tempSubjectDir)

def main():
  import argparse
  parser = argparse.ArgumentParser(
    usage="%(prog)s --subjects <LIDC_subjectID>\n\n"
    "This program will parse the DICOM and XML data for LIDC subject specified and generate"
    "DICOM representation for the segmentations and evaluations of the segmented nodule."
    "More details in a document to follow")
  parser.add_argument(
    '--subject-range',
    dest = "subjectRange",
    nargs=2,
    type=int,
    help = "Range of subject identifiers to be processed. Overrides individual subjects specified."
  )
  parser.add_argument(
    '--all-subjects',
    dest = "allSubjects",
    action="store_true",
    help = "Process all subjects (up to 1012). Overrides all other subject specifications."
  )
  parser.add_argument(
    '--subjects',
    type=int,
    nargs = '+',
    dest="subjectIDs",
    help='Identifier(s) (separated by space) of the subject to be processed.')
  parser.add_argument(
    '--log',
    dest="logFile",
    help="Location of the file to store processing log."
  )
  parser.add_argument(
    '--output-dir',
    dest="outputDir",
    help="Directory for storing the results of conversion."
  )

  args = parser.parse_args()

  if args.logFile:
    logging.basicConfig(filename=args.logFile,level=logging.INFO)
  else:
    logging.basicConfig(level=logging.INFO)
  logger = logging.getLogger("lidc2dicom")

  converter = LIDC2DICOMConverter()

  if args.outputDir:
    converter.tempDir = args.outputDir

  if args.subjectIDs:
    logger.info("Processing subjects "+str(args.subjectIDs))
    for s in args.subjectIDs:
      converter.convertForSubject(s)
  elif len(args.subjectRange):
    logger.info("Processing subjects from "+str(args.subjectRange[0])+" to "+str(args.subjectRange[1])+" inclusive")
    if args.subjectRange[1]<args.subjectRange[0]:
      logger.error("Invalid range.")
    for s in range(args.subjectRange[0],args.subjectRange[1]+1,1):
      converter.convertForSubject(s)
  elif args.allSubjects:
    logging.info("Processing all subjects from 1 to 1012.")
    for s in range(1,1013,1):
      converter.convertForSubject(s)


if __name__ == "__main__":
  main()

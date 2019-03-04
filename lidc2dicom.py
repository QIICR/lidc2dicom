from pathlib import Path
import lidc_conversion_utils.helpers as lidc_helpers
import os, itk, tempfile, json, pydicom, tempfile, shutil, sys
import subprocess
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
    self.tempDir="/Users/fedorov/Temp/LIDC_conversion2"

    # read GenericColors
    self.colors = []
    with open(self.colorsFile,'r') as f:
      for l in f:
        if l.startswith('#'):
          continue
        self.colors.append([int(c) for c in l.split(' ')[2:5]])

    self.conceptsDictionary = {}
    self.valuesDictionary = {}
    with open("concepts_dict.json") as cf:
      self.conceptsDictionary = json.load(cf)
    with open("values_dict.json") as vf:
      self.valuesDictionary = json.load(vf)

  def cleanUpTempDir(self, dir):
    shutil.rmtree(os.path.join(dir,"dicom2nrrd"))
    for p in Path(dir).glob("*.nrrd"):
      p.unlink()

  def convertSingleAnnotation(self, nCount, aCount, a, ctDCM, noduleUID, volume, seriesDir):
    with open(self.segTemplate,'r') as f:
      segJSON = json.load(f)

    # update as necessary!
    noduleName = "Nodule "+str(nCount+1)
    segName = "Nodule "+str(nCount+1) +" - Annotation " + a._nodule_id
    segJSON["segmentAttributes"][0][0]["SegmentDescription"] = segName
    segJSON["segmentAttributes"][0][0]["SegmentLabel"] = segName
    segJSON["SeriesDescription"] = "Segmentation of "+segName

    self.instanceCount = self.instanceCount+1
    if ctDCM.SeriesNumber != '':
      segJSON["SeriesNumber"] = str(int(ctDCM.SeriesNumber)+self.instanceCount)
    else:
      segJSON["SeriesNumber"] = str(self.instanceCount)

    for ci in range(3):
        segJSON["segmentAttributes"][0][0]["recommendedDisplayRGBValue"][ci] = self.colors[aCount+1][ci]

    segJSON["segmentAttributes"][0][0]["TrackingIdentifier"] = noduleName
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

    converterCmd = ['itkimage2segimage', "--inputImageList", nrrdSegFile, "--inputDICOMDirectory", seriesDir, "--inputMetadata", jsonSegFile, "--outputDICOM", dcmSegFile]
    self.logger.info("Converting to DICOM SEG with "+str(converterCmd))

    sp = subprocess.Popen(converterCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = sp.communicate()
    self.logger.info("itkimage2segimage stdout: "+stdout.decode('ascii'))
    self.logger.warning("itkimage2segimage stderr: "+stderr.decode('ascii'))

    segUID = None
    ctSeriesUID = None
    try:
      segDcm = pydicom.read_file(dcmSegFile)
      segUID = segDcm.SOPInstanceUID
      ctSeriesUID = segDcm.ReferencedSeriesSequence[0].SeriesInstanceUID
    except:
      self.logger.error("Failed to read Segmentation file")
      return

    with open(self.srTemplate,'r') as f:
      srJSON = json.load(f)

    srName = segName+" evaluations"
    srJSON["SeriesDescription"] = srName

    # be explicit about reader being anonymous
    srJSON["observerContext"] = {}
    srJSON["observerContext"]["ObserverType"] = "PERSON"
    srJSON["observerContext"]["PersonObserverName"] = "anonymous"

    self.instanceCount = self.instanceCount+1
    if ctDCM.SeriesNumber != '':
      srJSON["SeriesNumber"] = str(int(ctDCM.SeriesNumber)+self.instanceCount)
    else:
      srJSON["SeriesNumber"] = str(self.instanceCount)

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

    for attribute in self.conceptsDictionary.keys():
      # print(attribute+': '+str(getattr(a, attribute)))
      try:
        qItem = {}
        qItem["conceptCode"] = self.conceptsDictionary[attribute]
        qItem["conceptValue"] = self.valuesDictionary[attribute][str(getattr(a, attribute))]
        qualitativeEvaluations.append(qItem)
      except KeyError:
        self.logger.info("Skipping invalid attribute: "+attribute+': '+str(getattr(a, attribute)))
        continue

    srJSON["Measurements"][0]["measurementItems"] = measurementItems
    srJSON["Measurements"][0]["qualitativeEvaluations"] = qualitativeEvaluations
    srJSON["Measurements"][0]["segmentationSOPInstanceUID"] = segUID
    srJSON["Measurements"][0]["SourceSeriesForImageSegmentation"] = ctSeriesUID

    srJSON["Measurements"][0]["TrackingIdentifier"] = noduleName
    srJSON["Measurements"][0]["TrackingUniqueIdentifier"] = noduleUID

    srName = "Nodule "+str(nCount+1) +" - Annotation " + a._nodule_id + " measurements"
    jsonSRFile = os.path.join(self.tempSubjectDir,srName+'.json')
    with open(jsonSRFile, "w") as f:
      json.dump(srJSON, f, indent=2)

    dcmSRFile = os.path.join(self.tempSubjectDir,srName+'.dcm')
    converterCmd = ['tid1500writer', "--inputMetadata", jsonSRFile, "--inputImageLibraryDirectory", seriesDir, "--inputCompositeContextDirectory", self.tempSubjectDir, "--outputDICOM", dcmSRFile]
    self.logger.info("Converting with "+str(converterCmd))

    sp = subprocess.Popen(converterCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = sp.communicate()
    self.logger.info("tid1500writer stdout: "+stdout.decode('ascii'))
    self.logger.warning("tid1500writer stderr: "+stderr.decode('ascii'))

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

        sp = subprocess.Popen(plastimatchCmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        (stdout, stderr) = sp.communicate()
        self.logger.info("plastimatch stdout: "+stdout.decode('ascii'))
        self.logger.warning("plastimatch stderr: "+stderr.decode('ascii'))

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
    root = logging.getLogger()
    logging.basicConfig(filename=args.logFile,level=logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)s: %(levelname)s: %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)
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
  elif args.subjectRange is not None and len(args.subjectRange):
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

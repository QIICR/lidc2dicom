import itk, os, tempfile
import pandas as pd
import numpy as np
from subprocess import call
import pydicom

ITK_COORDINATE_UNKNOWN = 0
ITK_COORDINATE_Right = 2
ITK_COORDINATE_Left = 3
ITK_COORDINATE_Posterior = 4
ITK_COORDINATE_Anterior = 5
ITK_COORDINATE_Inferior = 8
ITK_COORDINATE_Superior = 9

ITK_COORDINATE_PrimaryMinor = 0
ITK_COORDINATE_SecondaryMinor = 8
ITK_COORDINATE_TertiaryMinor = 16

ITK_COORDINATE_ORIENTATION_RAS = ( ITK_COORDINATE_Right \
                            << ITK_COORDINATE_PrimaryMinor ) \
                                + ( ITK_COORDINATE_Anterior << ITK_COORDINATE_SecondaryMinor ) \
                                + ( ITK_COORDINATE_Superior << ITK_COORDINATE_TertiaryMinor )

ITK_COORDINATE_ORIENTATION_LPS = ( ITK_COORDINATE_Left \
                              << ITK_COORDINATE_PrimaryMinor ) \
                                + ( ITK_COORDINATE_Posterior << ITK_COORDINATE_SecondaryMinor ) \
                                + ( ITK_COORDINATE_Superior << ITK_COORDINATE_TertiaryMinor )

ITK_COORDINATE_ORIENTATION_LSP = ( ITK_COORDINATE_Left \
                                  << ITK_COORDINATE_PrimaryMinor ) \
                                + ( ITK_COORDINATE_Superior << ITK_COORDINATE_SecondaryMinor ) \
                                + ( ITK_COORDINATE_Posterior << ITK_COORDINATE_TertiaryMinor )


ITK_COORDINATE_ORIENTATION_RSP = (  ITK_COORDINATE_Right \
                                          <<  ITK_COORDINATE_PrimaryMinor ) \
                                          + ( ITK_COORDINATE_Superior << ITK_COORDINATE_SecondaryMinor ) \
                                       + ( ITK_COORDINATE_Posterior << ITK_COORDINATE_TertiaryMinor )

ITK_COORDINATE_ORIENTATION_RSA = ( ITK_COORDINATE_Right \
                                           << ITK_COORDINATE_PrimaryMinor ) \
                                    + ( ITK_COORDINATE_Superior << ITK_COORDINATE_SecondaryMinor ) \
                                       + ( ITK_COORDINATE_Anterior << ITK_COORDINATE_TertiaryMinor )

class SeriesGeometryChecker:

  # initialize with a pandas data frame containing series attributes
  # Columns that are expected to be present:
  #  - ImagePositionPatient
  #  - ImageOrientationPatient
  #  - SOPInstanceUID
  def __init__(self, series_df, warnings=False):
    self.series_df = series_df
    self.epsilon = 0.01
    self.computedSliceThickness = np.NaN
    self.warnings=warnings
    return
  #
  # now for each series and subseries, sort the images
  # by position and check for consistency
  #

  # TODO: more consistency checks:
  # - is there gantry tilt?
  # - are the orientations the same for all slices?
  def geometryOK(self):
    #
    # use the first file to get the ImageOrientationPatient for the
    # series and calculate the scan direction (assumed to be perpendicular
    # to the acquisition plane)
    #

    validGeometry = True
    ref = {}
    positions = []
    orientations = []
    for tag in ["ImagePositionPatient", "ImageOrientationPatient"]:
      for value in self.series_df[tag]:
        if not value or value == "NA":
          return False
        if tag == "ImagePositionPatient":
          #positions.append(np.array([float(zz) for zz in value.split('/')]))
          positions.append(np.array(value))
        if tag == "ImageOrientationPatient":
          #orientations.append(np.array([float(zz) for zz in value.split('/')]))
          orientations.append(np.array(value))

    # get the geometry of the scan
    # with respect to an arbitrary slice
    refPosition = positions[0]
    refOrientation = orientations[0]

    x = refOrientation[:3]
    y = refOrientation[3:]

    scanAxis = np.cross(x,y)

    #
    # for each file in series, calculate the distance along
    # the scan axis, sort files by this
    #
    sortList = []
    for p,o in zip(positions,orientations):
      dist = np.dot(p-refPosition, scanAxis)
      sortList.append((self.series_df["SOPInstanceUID"],dist))

    sortList = sorted(sortList, key=lambda x: x[1])
    '''
    for d in sortList:
      print(d[1])
    '''

    # confirm equal spacing between slices
    dist0 = sortList[1][1]-sortList[0][1]
    for i in range(len(sortList))[1:]:
      distN = sortList[i][1]-sortList[i-1][1]
      #print(str(sortList[i]))
      #print(distN)
      distDiff = distN-dist0
      if distDiff > self.epsilon:
        if self.warnings:
          print("Images are not equally spaced. Difference of %g in spacing detected!" % (distDiff))
        return False

    self.computedSliceThickness = dist0
    return True

def reconstructCTVolume(srcDir, destFile):
  tempDir = tempfile.mkdtemp()
  call(['dicom2nifti', srcDir, tempDir])
  print('dicom2nifti completed')
  outputFiles = os.listdir(tempDir)
  if not len(outputFiles):
    return False
  tempFile = os.path.join(tempDir, outputFiles[0])
  try:
    sitk.ReadImage(tempFile)
  except:
      return False
  try:
    # reorient
    ImageType = itk.Image[itk.SS, 3]
    reader = itk.ImageFileReader[ImageType].New()
    reader.SetFileName(tempFile)
    reader.Update()

    '''
    image = reader.GetOutput()
    reorient = itk.OrientImageFilter[ImageType,ImageType].New()
    reorient.SetInput(image)
    reorient.SetDesiredCoordinateOrientation(ITK_COORDINATE_ORIENTATION_RSP)
    reorient.Update();
    reoriented = reorient.GetOutput()
    '''

    writer = itk.ImageFileWriter[ImageType].New()
    writer.SetInput(reoriented)
    writer.SetFileName(destFile)
    writer.Update()

    return True
  except:
    return False

  return True

def checkSeriesGeometry(seriesDir):
  ipp = []
  iop = []
  uid = []
  print("Series dir: "+seriesDir)
  for f in os.listdir(seriesDir):
    try:
      dcm = pydicom.read_file(os.path.join(seriesDir,f))
    except:
      continue
    ipp.append(dcm.ImagePositionPatient)
    iop.append(dcm.ImageOrientationPatient)
    uid.append(dcm.SOPInstanceUID)
  series_df = pd.DataFrame(data={"ImagePositionPatient":ipp,"ImageOrientationPatient":iop,"SOPInstanceUID":uid})

  checker =     SeriesGeometryChecker(series_df)
  geometryOK = checker.geometryOK()

  if geometryOK:
    print("Geometry checks ok. Computed slice thickness: %f. DICOM SliceThickness: %s" % (checker.computedSliceThickness, dcm.SliceThickness))
    return True
  else:
    print("Geometry checks failed!")
    return False

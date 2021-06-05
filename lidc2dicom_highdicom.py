import argparse
from pathlib import Path
import logging
from typing import List
import os
import json
import sys

import pylidc as pl
import numpy as np

from pydicom import Dataset
from pydicom.sr.codedict import codes
from pydicom.valuerep import PersonName

from highdicom.version import __version__ as highdicom_version
from highdicom import AlgorithmIdentificationSequence, UID
from highdicom.seg import (
    Segmentation,
    SegmentDescription,
    SegmentAlgorithmTypeValues,
    SegmentationTypeValues
)

from highdicom.sr import (
    AlgorithmIdentification,
    CodedConcept,
    Comprehensive3DSR,
    FindingSite,
    Measurement,
    MeasurementReport,
    ObserverContext,
    ObservationContext,
    PersonObserverIdentifyingAttributes,
    QualitativeEvaluation,
    ReferencedSegment,
    RelationshipTypeValues,
    SourceImageForMeasurement,
    SourceSeriesForSegmentation,
    TrackingIdentifier,
    VolumetricROIMeasurementsAndQualitativeEvaluations,
)

import lidc_conversion_utils.helpers as lidc_helpers


class LIDC2DICOMConverter:

    def __init__(self, args):
        self.logger = logging.getLogger("lidc2dicom")

        self.args = args
        self.output_dir = args.output_dir
        self.uid_output_names = args.uid_output_names

        self.colors_file = "GenericColors.txt"

        # read GenericColors
        self.colors = []
        with open(self.colors_file, 'r') as f:
            for l in f:
                if l.startswith('#'):
                    continue
                self.colors.append([int(c) for c in l.split(' ')[2:5]])

        self.concepts_dictionary = {}
        self.values_dictionary = {}
        with open("concepts_dict.json") as cf:
            self.concepts_dictionary = json.load(cf)
        with open("values_dict.json") as vf:
            self.values_dictionary = json.load(vf)

        self.series_count = 1000

    def get_segment_description(
        self,
        segment_number: int,
        nodule_name: str,
        seg_name: str,
        nodule_uid: str,
        display_color: List[int]
    ):
        """Get a segment description object describing a segment of a segmentation image.

        In this application, each segment represents the annotation of a single nodule by
        a single reader. This method includes all information common to all segments,
        including the fact that this is a manual segmentation of a nodule found in the lung,
        using standard coding schemes.

        Parameters
        ----------
        segment_number: int
            Number of the segment within the segmentation image.
        nodule_name: str
            Name for the nodule this segment represents.
        seg_name: str
            Description of the segment.
        nodule_uid: str
            Tracking unique identifier for the nodule this
        display_color: List[int]
            Suggested display color for this segment when rendered by viewers.

        Returns
        -------
        highdicom.sr.SegmentDescription
            Description of the segment.

        """
        # Descriptive information about this segment
        seg_desc = SegmentDescription(
            segment_number=segment_number,
            segment_label=seg_name,
            segmented_property_category=codes.SCT.MorphologicallyAbnormalStructure,
            segmented_property_type=codes.SCT.Nodule,
            algorithm_type=SegmentAlgorithmTypeValues.MANUAL,
            tracking_uid=nodule_uid,
            tracking_id=nodule_name,
            anatomic_regions=[codes.SCT.Lung],
        )
        seg_desc.SegmentDescription = seg_name
        seg_desc.RecommendedDisplayCIELabValue = display_color

        return seg_desc

    def get_segmentation_dataset(
        self,
        ct_datasets: List[Dataset],
        pixel_array: np.ndarray,
        seg_descs: List[SegmentDescription],
        series_number: int,
        series_description: str
    ):
        """Construct the SOP Instance dataset for a segmentation image.

        This method includes information common to all segmentation images,
        such as the manufacturer name and software version.

        Parameters
        ----------
        ct_datasets: List[pydicom.dataset.Dataset]
            List of CT datasets from which the segmentations were derived.
        pixel_array: np.ndarray
            Pixel array of the segmentation masks
        seg_descs: List[SegmentDescription]
            Descriptions of each segment in the segmentation image.
        series_number: int
            Series number to apply to the newly created segmentation image.
        series_description: str
            Series description of the newly-created segmentation image.

        Returns
        -------
        highdicom.seg.Segmentation
            Dataset for the newly-created Segmentation SOP Instance.

        """
        anonymous_person_name = PersonName.from_named_components(
            given_name='anonymous',
            family_name='observer'
        )
        seg_dcm = Segmentation(
            source_images=ct_datasets,
            pixel_array=pixel_array,
            segmentation_type=SegmentationTypeValues.BINARY,
            segment_descriptions=seg_descs,
            series_instance_uid=UID(),
            series_number=series_number,
            sop_instance_uid=UID(),
            instance_number=1,
            manufacturer="highdicom developers",
            manufacturer_model_name="highdicom",
            software_versions=f"{highdicom_version}",
            device_serial_number='1',
            content_description="Lung nodule segmentation",
            content_creator_name=anonymous_person_name,
            series_description=series_description
        )

        # Add in some extra information
        seg_dcm.BodyPartExamined = "LUNG"
        seg_dcm.ClinicalTrialSeriesID = "Session1"
        seg_dcm.ClinicalTrialTimePointID = "1"
        seg_dcm.ClinicalTrialCoordinatingCenterName = "TCIA"
        seg_dcm.ContentLabel = "SEGMENTATION"

        return seg_dcm

    def get_roi_measurements_and_evaluations(
        self,
        ann: pl.Annotation,
        ct_datasets: List[Dataset],
        seg_dcm: Dataset,
        segment_number: int,
        nodule_uid: str,
        nodule_name: str
    ):
        """Construct a measurements group for a given annotation for inclusion
        in the SR document content tree.

        In this application, a measurements group is used to encode
        measurements derived from a single segment in the segmentation image,
        i.e. an annotation of a single nodule by a single reader. The
        measurements consist of the volume, surface area, and diameter of the
        nodule, as well as several qualitative evaluations that record rating
        of the subtlety, internal structure, calcification, sphericity, margin,
        lobulation, spiculation, texture, and malignancy of the nodule. This
        information is provided directly from the pylidc.Annotation object.

        Parameters
        ----------
        ann: pl.Annotation
            The annotation object, in which all measurements and evaluations
            are recorded.
        ct_datasets: List[pydicom.dataset.Dataset]
            List of CT datasets from which the segmentations were derived.
        seg_dcm: highdicom.seg.Segmentation
            Segmentation image dataset containing the segmentations that will
            be referenced in the SR content.
        segment_number: int
            Number of the segment within the segmentation image that
            corresponds to this annotation.
        nodule_uid: str
            Tracking unique identifier for the nodule this
        nodule_name: str
            Name for the nodule this segment represents.

        Returns
        -------
        highdicom.sr.VolumetricROIMeasurementsAndQualitativeEvaluations
            SR template (TID 1411) object representing measurements and
            evaluations for a single nodule as annotated by a single reader.

        """
        # Get measurements from a single annotation and encode in a TID 1411
        # template

        # Identify pylidc as the "algorithm" creating the annotations
        pylidc_algo_id = AlgorithmIdentification(name='pylidc', version=pl.__version__)

        # Describe the anatomic site at which observations were made
        finding_sites = [FindingSite(anatomic_location=codes.SCT.Lung)]

        referenced_segment = ReferencedSegment(
            sop_class_uid=seg_dcm.SOPClassUID,
            sop_instance_uid=seg_dcm.SOPInstanceUID,
            segment_number=segment_number,
            source_series=SourceSeriesForSegmentation.from_source_image(
                ct_datasets[0]
            )
        )

        # Describe the imaging measurements for the image region defined above
        referenced_images = [
            SourceImageForMeasurement.from_source_image(ds)
            for ds in ct_datasets
        ]
        # Volume measurement
        volume_measurement = Measurement(
            name=codes.SCT.Volume,
            tracking_identifier=TrackingIdentifier(uid=UID()),
            value=ann.volume,
            unit=codes.UCUM.CubicMillimeter,
            referenced_images=referenced_images,
            algorithm_id=pylidc_algo_id
        )
        # Diameter measurement
        diameter_measurement = Measurement(
            name=codes.SCT.Diameter,
            tracking_identifier=TrackingIdentifier(uid=UID()),
            value=ann.diameter,
            unit=codes.UCUM.Millimeter,
            referenced_images=referenced_images,
            algorithm_id=pylidc_algo_id
        )
        # Surface area measurement
        surface_area_measurement = Measurement(
            name=CodedConcept(value='C0JK', scheme_designator='IBSI', meaning="Surface area of mesh"),
            tracking_identifier=TrackingIdentifier(uid=UID()),
            value=ann.surface_area,
            unit=codes.UCUM.SquareMillimeter,
            referenced_images=referenced_images,
            algorithm_id=pylidc_algo_id
        )

        # Qualitative evaluations
        qualitative_evaluations = []
        for attribute in self.concepts_dictionary.keys():
            try:
                qualitative_evaluations.append(
                    QualitativeEvaluation(
                        name=CodedConcept(**self.concepts_dictionary[attribute]),
                        value=CodedConcept(**self.values_dictionary[attribute][str(getattr(ann, attribute))])
                    )
                )
            except KeyError:
                self.logger.info(f"Skipping invalid attribute: {attribute} {getattr(ann, attribute)}")
                continue

        # Compile into TID1411
        roi_measurements = VolumetricROIMeasurementsAndQualitativeEvaluations(
            tracking_identifier=TrackingIdentifier(
                uid=nodule_uid,
                identifier=nodule_name
            ),
            referenced_segment=referenced_segment,
            finding_type=codes.SCT.Nodule,
            measurements=[volume_measurement, diameter_measurement, surface_area_measurement],
            qualitative_evaluations=qualitative_evaluations,
            finding_sites=finding_sites
        )

        return roi_measurements

    def get_sr_dataset(
        self,
        roi_measurements: List[VolumetricROIMeasurementsAndQualitativeEvaluations],
        ct_datasets: List[Dataset],
        seg_dataset: Segmentation,
        series_number: int,
        series_description: str
    ):
        """Construct a Structured Report containing measurements and
        evaluations for one or more annotations.

        Depending on the command-line parameters, each SR document may contain
        measurements and evaluations for the annotation of a single nodule by a
        single reader, or in the case of "--composite", measurements and
        evaluations of one or more nodules, each annotated by one or more
        distinct readers.

        Parameters
        ----------
        roi_measurements: List[highdicom.sr.VolumetricROIMeasurementsAndQualitativeEvaluations]
            List of measurement groups, each containing to the measurements and
            evaluations of a single nodule by a single reader, and referencing
            a single segment in a segmentation image.
        ct_datasets: List[pydicom.dataset.Dataset]
            List of CT datasets from which the segmentations were derived.
        seg_dataset: Segmentation
            Segmentation image referenced by the SR.
        series_number: int
            Series number of the newly created SR document.
        series_description: str
            Series description of the newly-created SR document.

        Returns
        -------
        highdicom.sr.Comprehensive3DSR
            Structured Report document dataset using template TID 1500
            containing measurements and evaluations for one or more nodules.

        """
        # Be explicit about reader being anonymous
        anonymous_person_name = PersonName.from_named_components(
            given_name='anonymous',
            family_name='observer'
        )
        observer_context = ObserverContext(
            observer_type=codes.DCM.Person,
            observer_identifying_attributes=PersonObserverIdentifyingAttributes(
                name=anonymous_person_name
            )
        )
        observation_context = ObservationContext(
            observer_person_context=observer_context
        )

        measurement_report = MeasurementReport(
            observation_context=observation_context,
            procedure_reported=codes.LN.CTUnspecifiedBodyRegion,
            imaging_measurements=roi_measurements
        )

        # Create the Structured Report instance
        sr_dcm = Comprehensive3DSR(
            evidence=ct_datasets + [seg_dataset],
            content=measurement_report[0],
            series_number=series_number,
            series_instance_uid=UID(),
            sop_instance_uid=UID(),
            instance_number=1,
            manufacturer='highdicom developers',
            is_complete=True,
            is_verified=True,
            verifying_observer_name=anonymous_person_name,
            verifying_organization='anonymous',
            series_description=series_description
        )

        return sr_dcm

    def convert_single_annotation(
        self,
        n_count: int,
        a_count: int,
        ann: pl.Annotation,
        ct_datasets: List[Dataset],
        nodule_uid: str,
        series_dir: str,
        scan: pl.Scan
    ):
        """Convert a single annotation into a DICOM Segmentation image and a
        DICOM SR document.

        A single annotation is a segmentation and a collection of related
        measurements and evaluations created by a single reader for a single
        nodule. The resulting DICOM datasets are saved to the filesystem as
        files.

        Parameters
        ----------
        n_count: int
            Nodule count, i.e. running count of nodules within this scan.
            This is used to name the nodule.
        a_count: int
            Annotation count, i.e. running count of annotations for this
            nodule. This is used to name the annotation.
        ann: pylidc.Annotation
            Pylidc Annotation object containing all the information for this
            annotation, including the segmentation itself as well as
            corresponding measurements and evaluations.
        ct_datasets: List[pydicom.dataset.Dataset]
            List of CT datasets from which the annotations were derived.
        nodule_uid: str
            Tracking unique identifier for the nodule this.
        series_dir: str
            Path to directory containing the files for the original CT DICOM
            series.
        scan: pylidc.Scan
            Pylidc Scan object containing all information about the CT scan
            on which the annotation was performed.

        """
        nodule_name = f"Nodule {n_count + 1}"
        seg_name = f"Nodule {n_count + 1} - Annotation {ann._nodule_id}"

        # Choose series numbers for the new series and increment the counter
        seg_series_number = self.series_count
        sr_series_number = self.series_count + 1
        self.series_count += 2

        self.logger.info("Creating DICOM SEG")

        # Construct an empty mask the same size as the input series
        image_size = (ct_datasets[0].Rows, ct_datasets[0].Columns, len(ct_datasets))
        mask = np.zeros(image_size, np.uint8)

        # Fill in the mask elements with the segmentation
        mask[ann.bbox()] = ann.boolean_mask().astype(np.int8)

        # Find the subset of the source images relevant for the segmentation
        ct_subset = ct_datasets[ann.bbox()[2]]
        mask_subset = mask[(slice(None), slice(None), ann.bbox()[2])]
        mask_subset = np.moveaxis(mask_subset, 2, 0)

        seg_desc = self.get_segment_description(
            segment_number=1,
            nodule_name=nodule_name,
            seg_name=seg_name,
            nodule_uid=nodule_uid,
            display_color=self.colors[a_count + 1]
        )

        seg_dcm = self.get_segmentation_dataset(
            ct_datasets=ct_subset,
            pixel_array=mask_subset,
            seg_descs=[seg_desc],
            series_number=seg_series_number,
            series_description=f"Segmentation of {seg_name}"
        )

        # Save the file
        seg_series_dir = os.path.join(self.subject_dir, seg_dcm.SeriesInstanceUID)
        os.makedirs(seg_series_dir)
        if self.uid_output_names:
            dcm_seg_file = os.path.join(seg_series_dir, seg_dcm.SOPInstanceUID + '.dcm')
        else:
            dcm_seg_file = os.path.join(seg_series_dir, seg_name + '.dcm')
        seg_dcm.save_as(dcm_seg_file)

        self.logger.info("Creating DICOM SR")

        sr_name = f"Nodule {n_count + 1} - Annotation {ann._nodule_id} measurements"

        roi_measurements = self.get_roi_measurements_and_evaluations(
            ann=ann,
            ct_datasets=ct_subset,
            seg_dcm=seg_dcm,
            segment_number=1,
            nodule_uid=nodule_uid,
            nodule_name=nodule_name
        )

        sr_dcm = self.get_sr_dataset(
            roi_measurements=[roi_measurements],
            ct_datasets=ct_subset,
            seg_dataset=seg_dcm,
            series_number=sr_series_number,
            series_description=sr_name
        )

        # Save the file
        sr_series_dir = os.path.join(self.subject_dir, sr_dcm.SeriesInstanceUID)
        os.makedirs(sr_series_dir)
        if self.uid_output_names:
            dcm_sr_file = os.path.join(sr_series_dir, sr_dcm.SOPInstanceUID + '.dcm')
        else:
            dcm_sr_file = os.path.join(sr_series_dir, sr_name + '.dcm')
        sr_dcm.save_as(dcm_sr_file)

    def convert_for_scan(
        self,
        scan: pl.Scan, 
        ct_datasets: List[Dataset],
        series_dir: str
    ):
        """Convert all annotations within a given scan to DICOM format as
        single Segmentations and SRs.

        Resulting DICOM files are stored to the filesystem.

        Parameters
        ----------
        scan: pylidc.Scan
            Pylidc Scan object containing all information about the CT scan
            on which the annotation was performed.
        ct_datasets: List[pydicom.dataset.Dataset]
            List of CT datasets from which the annotations were derived.
        series_dir: str
            Path to directory containing the files for the original CT DICOM
            series.

        """
        # Iterate over all nodules available for this subject
        anns = scan.annotations
        self.logger.info(f'Have {len(anns)} annotations for subject {scan.patient_id}')

        self.instance_count = 0

        clustered_annotation_ids = []

        for n_count, nodule in enumerate(scan.cluster_annotations()):
            nodule_uid = UID() 

            for a_count, ann in enumerate(nodule):
                clustered_annotation_ids.append(ann.id)
                self.convert_single_annotation(
                    n_count=n_count,
                    a_count=a_count,
                    ann=ann,
                    ct_datasets=ct_datasets,
                    nodule_uid=nodule_uid,
                    series_dir=series_dir,
                    scan=scan
                )

        if len(clustered_annotation_ids) != len(anns):
            n_missing = len(anns) - len(clustered_annotation_ids)
            self.logger.warning(
                f"{n_missing} annotations unaccounted for!"
            )

        for ua in anns:
            if ua.id not in clustered_annotation_ids:
                a_count = a_count + 1
                n_count = n_count + 1
                nodule_uid = UID()
                self.convert_single_annotation(
                    n_count=n_count,
                    a_count=a_count,
                    ann=ua,
                    ct_datasets=ct_datasets,
                    nodule_uid=nodule_uid,
                    series_dir=series_dir,
                    scan=scan
                )

    def convert_for_scan_composite(
        self,
        scan: pl.Scan,
        ct_datasets: List[Dataset],
        series_dir: str
    ):
        """Convert all annotations for a scan into a single "composite" DICOM
        Segmentation image and accompanying DICOM SR document.

        A single annotation is a segmentation and a collection of related
        measurements and evaluations created by a single reader for a single
        nodule. The resulting DICOM datasets are saved to the filesystem as
        files. This method combines all such annotations for a single scan
        into a single DICOM Segmentation image and SR document (in contrast
        to the convert_for_scan method).

        Parameters
        ----------
        scan: pylidc.Scan
            Pylidc Scan object containing all information about the CT scan
            on which the annotation was performed.
        ct_datasets: List[pydicom.dataset.Dataset]
            List of CT datasets from which the annotations were derived.
        series_dir: str
            Path to directory containing the files for the original CT DICOM
            series.

        """
        n_annotations = len(scan.annotations)
        self.logger.info(
            f'Have {n_annotations} annotations for subject {scan.patient_id}'
        )
        if n_annotations == 0:
            # Nothing to do
            return

        # Choose series numbers for the new series and increment the counter
        seg_series_number = self.series_count
        sr_series_number = self.series_count + 1
        self.series_count += 2

        image_size = (
            ct_datasets[0].Rows,
            ct_datasets[0].Columns,
            len(ct_datasets)
        )

        total_ann_ind = 0

        all_roi_measurements = []
        for n_count, nodule in enumerate(scan.cluster_annotations()):
            nodule_uid = UID()
            nodule_name = f"Nodule {n_count + 1}"

            for a_count, a in enumerate(nodule):

                seg_name = f"Nodule {n_count + 1} - Annotation {a._nodule_id}"

                # A number for this segment within the segmentation object
                # (must begin at 1)
                segment_number = total_ann_ind + 1

                seg_desc = self.get_segment_description(
                    segment_number=segment_number,
                    nodule_name=nodule_name,
                    seg_name=seg_name,
                    nodule_uid=nodule_uid,
                    display_color=self.colors[total_ann_ind + 1]
                )

                # Construct an empty mask the same size as the input series
                mask = np.zeros(image_size, np.uint8)

                # Fill in the mask elements with the segmentation
                mask[a.bbox()] = a.boolean_mask().astype(np.int8)
                mask = np.moveaxis(mask, 2, 0)

                if total_ann_ind == 0:
                    # Need to create the segmentation object
                    seg_dcm = self.get_segmentation_dataset(
                        ct_datasets=ct_datasets,
                        pixel_array=mask,
                        seg_descs=[seg_desc],
                        series_number=seg_series_number,
                        series_description='Segmentation of All Nodules'
                    )
                else:
                    # Add new segment to the existing object
                    seg_dcm.add_segments(mask, [seg_desc])

                roi_measurements = self.get_roi_measurements_and_evaluations(
                    ann=a,
                    ct_datasets=ct_datasets,
                    seg_dcm=seg_dcm,
                    segment_number=segment_number,
                    nodule_uid=nodule_uid,
                    nodule_name=nodule_name
                )
                all_roi_measurements.append(roi_measurements)

                total_ann_ind += 1

        # Save the file
        seg_series_dir = os.path.join(self.subject_dir, seg_dcm.SeriesInstanceUID)
        os.makedirs(seg_series_dir)
        if self.uid_output_names:
            dcm_seg_file = os.path.join(seg_series_dir, seg_dcm.SOPInstanceUID + '.dcm')
        else:
            dcm_seg_file = os.path.join(seg_series_dir, 'all_segmentations.dcm')
        seg_dcm.save_as(dcm_seg_file)

        sr_dcm = self.get_sr_dataset(
            roi_measurements=all_roi_measurements,
            ct_datasets=ct_datasets,
            seg_dataset=seg_dcm,
            series_number=sr_series_number,
            series_description='All nodules measurements'
        )

        # Save the file
        sr_series_dir = os.path.join(self.subject_dir, sr_dcm.SeriesInstanceUID)
        os.makedirs(sr_series_dir)
        if self.uid_output_names:
            dcm_sr_file = os.path.join(sr_series_dir, sr_dcm.SOPInstanceUID + '.dcm')
        else:
            dcm_sr_file = os.path.join(sr_series_dir, 'all_measurements.dcm')
        sr_dcm.save_as(dcm_sr_file)

    def convert_for_subject(self, subject_id: int, composite: bool = False):
        """Convert all scans for a given subject to DICOM format.

        Resulting DICOM files are stored to the filesystem.

        Parameters
        ----------
        subject_id: int
            LIDC subject ID (between 1 and 1012 inclusive).
        composite: bool
            If True, create Segmentation and SR objects containing all
            annotations for each scan combined. If False (the default),
            a single SR and Segmentation are created for each annotation.

        """
        s = 'LIDC-IDRI-%04i' % subject_id
        self.logger.info("Processing subject %s" % (s))
        scans = pl.query(pl.Scan).filter(pl.Scan.patient_id == s)
        self.logger.info(" Found %d scans" % (scans.count()))

        for scan in scans:
            study_uid = scan.study_instance_uid
            series_uid = scan.series_instance_uid
            series_dir = Path(scan.get_path_to_dicom_files())
            if not os.path.exists(series_dir):
                self.logger.error("Files not found for subject " + s)
                return

            # Reset the series counter (used to number new series)
            self.series_count = 100

            try:
                ct_datasets = scan.load_all_dicom_images()
            except Exception:
                self.logger.error("Failed to read input CT files")
                return

            ok = lidc_helpers.checkSeriesGeometry(str(series_dir))
            if not ok:
                self.logger.warning("Geometry inconsistent for subject %s" % (s))

            self.subject_dir = os.path.join(self.output_dir, s, study_uid)
            os.makedirs(self.subject_dir, exist_ok=True)

            if composite:
                self.convert_for_scan_composite(scan, ct_datasets, series_dir)
            else:
                self.convert_for_scan(scan, ct_datasets, series_dir)


def main():
    """Main entrypoint function. Parses command-line arguments, creates converter object
    and initiates conversion.

    """
    parser = argparse.ArgumentParser(
        usage="%(prog)s --subjects <LIDC_subjectID>\n\n"
        "This program will parse the DICOM and XML data for LIDC subject specified and generate"
        "DICOM representation for the segmentations and evaluations of the segmented nodule."
        "More details in a document to follow"
    )
    parser.add_argument(
        '--subject-range',
        '-r',
        dest="subject_range",
        nargs=2,
        type=int,
        help="Range of subject identifiers to be processed. Overrides individual subjects specified."
    )
    parser.add_argument(
        '--all-subjects',
        '-a',
        dest="all_subjects",
        action="store_true",
        help="Process all subjects (up to 1012). Overrides all other subject specifications."
    )
    parser.add_argument(
        '--subjects',
        '-s',
        type=int,
        nargs='+',
        dest="subject_ids",
        help='Identifier(s) (separated by space) of the subject to be processed.'
    )
    parser.add_argument(
        '--log',
        '-l',
        dest="log_file",
        help="Location of the file to store processing log."
    )
    parser.add_argument(
        '--output-dir',
        '-o',
        dest="output_dir",
        help="Directory for storing the results of conversion."
    )
    parser.add_argument(
        '--composite',
        '-c',
        action="store_true",
        default=False,
        dest="composite",
        help="Make composite objects (1 SEG and 1 SR that contain all segmentations/measurement for "
             "all nodes/annotations). Composite objects will not be generated by default."
    )
    parser.add_argument(
        '--images-dir',
        '-i',
        dest="images_dir",
        help="Directory with the CT images of the LIDC-IDRI collection. The directory should be organized "
             "following this pattern: <subject ID>/<study UID>/<series UID>."
    )
    parser.add_argument(
        '--uid-output-names',
        '-u',
        action='store_true',
        dest="uid_output_names",
        help="Name output DICOM files by their SOP Instance UID. If False, a human-readable name is used."
    )

    args = parser.parse_args()

    if args.log_file:
        root = logging.getLogger()
        logging.basicConfig(filename=args.log_file, level=logging.INFO)
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(name)s: %(levelname)s: %(message)s')
        handler.setFormatter(formatter)
        root.addHandler(handler)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("lidc2dicom")

    converter = LIDC2DICOMConverter(args)

    if args.output_dir:
        converter.output_dir = args.output_dir

    if args.subject_ids:
        logger.info(f"Processing subjects {args.subject_ids}")
        for s in args.subject_ids:
            converter.convert_for_subject(s, composite=args.composite)
    elif args.subject_range is not None and len(args.subject_range):
        logger.info(f"Processing subjects from {args.subject_range[0]} to {args.subject_range[1]} inclusive")
        if args.subject_range[1] < args.subject_range[0]:
            logger.error("Invalid range.")
        for s in range(args.subject_range[0], args.subject_range[1] + 1, 1):
            converter.convert_for_subject(s, composite=args.composite)
    elif args.all_subjects:
        logging.info("Processing all subjects from 1 to 1012.")
        for s in range(1, 1013, 1):
            converter.convert_for_subject(s, composite=args.composite)


if __name__ == "__main__":
    main()

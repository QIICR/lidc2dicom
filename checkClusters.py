import lidc_conversion_utils.helpers as lidc_helpers
import os, itk, tempfile, json, pydicom, tempfile, shutil
from subprocess import call
import pylidc as pl
import numpy as np
import glob
import logging
import pylidc as pl

def checkSubject(subjectID, show):
  s = 'LIDC-IDRI-%04i' % subjectID
  print("Checking subject "+s)
  scan = pl.query(pl.Scan).filter(pl.Scan.patient_id == s).first()
  if scan is None:
    print("%s not found!" % (s))
    return
  nodules = scan.cluster_annotations()
  annotations = pl.query(pl.Annotation).join(pl.Scan).filter(pl.Scan.patient_id == s)
  print("%i nodules and %d annotations" % (len(nodules), annotations.count()))
  annotationsInNodules = 0
  annotationsInNodulesList = []
  annotationsNotInNodulesList = []
  for nCount,nodule in enumerate(nodules):
    print("  Nodule %d has %d annotations" % (nCount+1, len(nodule)))
    for a in nodule:
      annotationsInNodulesList.append(a.id)
  if len(annotationsInNodulesList) != annotations.count():
    print("   WARNING: %d annotations unaccounted for!" % (annotations.count()-len(annotationsInNodulesList)))
    allAnnotationsList = []
    for a in annotations:
      if a.id not in annotationsInNodulesList:
        annotationsNotInNodulesList.append(a.id)
        print("%d (%s) not assigned to a nodule" % (a.id, a._nodule_id))
      else:
        print("%d (%s) assigned to a nodule" % (a.id, a._nodule_id))
  annotations2show = []

  if show == 'all':
    annotations2show = [[a] for a in annotations]
  elif show == 'missing':
    annotations2show = [[a] for a in annotations if a.id in annotationsNotInNodulesList]
  elif show != None:
    annotations2show = [[a] for a in annotations if a._nodule_id == show]

  if show != None and len(annotations2show) == 0:
    print("OOPS: Could not figure out what to show with the parameters given.")
  if len(annotations2show) != 0:
    scan.visualize(annotation_groups=annotations2show)

def main():
  import argparse
  parser = argparse.ArgumentParser(
    usage="%(prog)s --subjects <LIDC_subjectID>\n\n"
    "This program will check status of annotations and node clusters.")
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
    '--show',
    type=str,
    dest="showAnnotations",
    help="Show annotations. Allowed choices \'all\', \'missing\', or _nodule_id of the specific annotation to show."
  )

  args = parser.parse_args()

  if args.subjectIDs:
    logging.info("Processing subjects "+str(args.subjectIDs))
    for s in args.subjectIDs:
      checkSubject(s, show=args.showAnnotations)
  elif len(args.subjectRange):
    logging.info("Processing subjects from "+str(args.subjectRange[0])+" to "+str(args.subjectRange[1])+" inclusive")
    if args.subjectRange[1]<args.subjectRange[0]:
      print("Invalid range.")
    for s in range(args.subjectRange[0],args.subjectRange[1]+1,1):
      checkSubject(s, show=args.showAnnotations)
  elif args.allSubjects:
    logging.info("Processing all subjects from 1 to 1012.")
    for s in range(1,1013,1):
      checkSubject(s, show=args.showAnnotations)


if __name__ == "__main__":
  main()

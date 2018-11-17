import pylidc as pl


for subjectID in range(1,1013):
  s = 'LIDC-IDRI-%04i' % subjectID
  scans = pl.query(pl.Scan).filter(pl.Scan.patient_id == s)
  if scans.count()>1:
    print("%s has %d scans" % (s, scans.count()))
    for i,scan in enumerate(scans):
      print("  Scan %d has %d annotations" % (i+1,len(scan.annotations)))

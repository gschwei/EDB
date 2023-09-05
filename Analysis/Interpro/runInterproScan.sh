# created by Gabriel on 03/03/16

# Here, I predict InterproScan domains for the secretome of each species.

echo "run Interpro 1"
java -jar interproscan-5.jar -u . -appl ProSitePatterns,SUPERFAMILY,Hamap,Pfam,ProSiteProfiles,ProDom,TIGRFAM,PRINTS,SMART,PIRSF,PANTHER,SFLD,CDD,MobiDBLite -i ../"MergedSecretome.fa" -f tsv -o ../"InterproScanAll.txt"

echo "run Interpro 2"
java -jar interproscan-5.jar -u . -appl Gene3D,Coils -i ../"MergedSecretome.fa" -f tsv -o ../"InterproScan2All2.txt"

echo "Interpro is done. Bye."

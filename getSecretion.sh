# created by Gabriel on 10/02/17
# last modified by Gabriel on 10/02/17

# Here, I predict secretion with SignalP 4.1, TMHMM 2.0c and store the results.
# NB: there is a problem with running phobius; therefore, I use this program in its online version.

for g in "Mepe" "Usho" "Usbr" "Usma" "Spsc" "Sprz" "Sprs"
do   
    mkdir $g
    #echo $g': run Phobius'
    #phobius -short ../../Data/$g/$g"_valid_prot.fa" > $g/$g"_phobius.txt"
    echo $g': run SignalP'
    signalp4.1 -f short ../../Data/$g/$g"_valid_prot.fa" > $g/$g"_signalp.txt" 
    echo $g': run TMHMM'
    /home/schweizerg/Programs/tmhmm/2.0c/bin/tmhmm -short ../../Data/$g/$g"_valid_prot.fa" > $g/$g"_TMHMM.txt"
    echo $g 'is done. Go to next species. Bye.'
done

rm Vertebrate*.fa
# Chicken
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Gallus_gallus:/" >> Vertebrate1.fa
# Cod
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/gadus_morhua/dna/Gadus_morhua.gadMor1.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Gadus_morhua:/" >> Vertebrate1.fa
# Cow
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Bos_taurus:/" >> Vertebrate1.fa
# Dog
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Canis_familiaris:/" >> Vertebrate2.fa
# Duck
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/anas_platyrhynchos/dna/Anas_platyrhynchos.BGI_duck_1.0.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Anas_platyrhynchos:/" >> Vertebrate2.fa
# Frog
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.JGI_4.2.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Xenopus_tropicalis:/" >> Vertebrate2.fa
# Horse
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/equus_caballus/dna/Equus_caballus.EquCab2.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Equus_caballus:/" >> Vertebrate3.fa
# Rabbit
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Oryctolagus_cuniculus:/" >> Vertebrate3.fa
# Sheep
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/ovis_aries/dna/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Ovis_aries:/" >> Vertebrate3.fa
# Wild boar (pig)
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa10.2.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Sus_scrofa:/" >> Vertebrate4.fa
# Tilapia
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/oreochromis_niloticus/dna/Oreochromis_niloticus.Orenil1.0.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Oreochromis_niloticus:/" >> Vertebrate4.fa
# Turkey
wget -qO- ftp://ftp.ensembl.org/pub/release-78/fasta/meleagris_gallopavo/dna/Meleagris_gallopavo.UMD2.dna_sm.toplevel.fa.gz | gunzip -c | sed "s/>/>Vertebrate:Meleagris_gallopavo:/" >> Vertebrate4.fa


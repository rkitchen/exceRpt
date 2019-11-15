cd /gpfs/scratch/fas/gerstein/rrk24/ANNOTATIONS/Genomes_BacteriaFungiMammalPlantProtistVirus/PlantGenomes_fasta

rm Plant*.fa

## Extract
gunzip *.fa.gz

## Fix missing chromosome number in first sequence of Solanum_lycopersicum.SL2.40.24.dna_sm.toplevel.fa
cat Solanum_lycopersicum.SL2.40.24.dna_sm.toplevel.fa | sed 's/^> dna/>0 dna/' > Solanum_lycopersicum.SL2.40.24.dna_sm.toplevel.fixed.fa 



cat Aegilops_tauschii.GCA_000347335.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Aegilops_tauschii:/" >> Plant1.fa
cat Amborella_trichopoda.GCA_000471905.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Amborella_trichopoda:/" >> Plant1.fa
cat Arabidopsis_lyrata.v.1.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Arabidopsis_lyrata:/" >> Plant1.fa
cat Arabidopsis_thaliana.TAIR10.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Arabidopsis_thaliana:/" >> Plant1.fa
cat Brachypodium_distachyon.v1.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Brachypodium_distachyon:/" >> Plant1.fa
cat Brassica_oleracea.v2.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Brassica_oleracea:/" >> Plant1.fa
cat Brassica_rapa.IVFCAASv1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Brassica_rapa:/" >> Plant1.fa

cat Chlamydomonas_reinhardtii.v3.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Chlamydomonas_reinhardtii:/" >> Plant2.fa
cat Cyanidioschyzon_merolae.ASM9120v1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Cyanidioschyzon_merolae:/" >> Plant2.fa
cat Glycine_max.V1.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Glycine_max:/" >> Plant2.fa
cat Hordeum_vulgare.082214v1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Hordeum_vulgare:/" >> Plant2.fa

cat Leersia_perrieri.Lperr_V1.4.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Leersia_perrieri:/" >> Plant3.fa
cat Medicago_truncatula.GCA_000219495.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Medicago_truncatula:/" >> Plant3.fa
cat Musa_acuminata.MA1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Musa_acuminata:/" >> Plant3.fa
cat Oryza_barthii.ABRL00000000.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_barthii:/" >> Plant3.fa
cat Oryza_brachyantha.Oryza_brachyantha.v1.4b.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_brachyantha:/" >> Plant3.fa
cat Oryza_glaberrima.AGI1.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_glaberrima:/" >> Plant3.fa
cat Oryza_glumaepatula.ALNU02000000.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_glumaepatula:/" >> Plant3.fa
cat Oryza_indica.ASM465v1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_indica:/" >> Plant3.fa
cat Oryza_meridionalis.ALNW00000000.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_meridionalis:/" >> Plant3.fa
cat Oryza_nivara.AWHD00000000.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_nivara:/" >> Plant3.fa
cat Oryza_punctata.AVCL00000000.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_punctata:/" >> Plant3.fa
cat Oryza_rufipogon.PRJEB4137.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_rufipogon:/" >> Plant3.fa
cat Oryza_sativa.IRGSP-1.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Oryza_sativa:/" >> Plant3.fa
cat Ostreococcus_lucimarinus.GCA_000092065.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Ostreococcus_lucimarinus:/" >> Plant3.fa
cat Physcomitrella_patens.ASM242v1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Physcomitrella_patens:/" >> Plant3.fa
cat Populus_trichocarpa.JGI2.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Populus_trichocarpa:/" >> Plant3.fa

cat Prunus_persica.GCA_000346465.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Prunus_persica:/" >> Plant4.fa
cat Selaginella_moellendorffii.v1.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Selaginella_moellendorffii:/" >> Plant4.fa
cat Setaria_italica.JGIv2.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Setaria_italica:/" >> Plant4.fa
cat Solanum_lycopersicum.SL2.40.24.dna_sm.toplevel.fixed.fa | sed "s/>/>Plant:Solanum_lycopersicum:/" >> Plant4.fa
cat Solanum_tuberosum.3.0.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Solanum_tuberosum:/" >> Plant4.fa
cat Sorghum_bicolor.Sorbi1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Sorghum_bicolor:/" >> Plant4.fa
cat Theobroma_cacao.GCA_000403535.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Theobroma_cacao:/" >> Plant4.fa
cat Vitis_vinifera.IGGP_12x.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Vitis_vinifera:/" >> Plant4.fa
cat Zea_mays.AGPv3.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Zea_mays:/" >> Plant4.fa

cat Triticum_aestivum.IWGSP1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Triticum_aestivum:/" >> Plant5.fa
cat Triticum_urartu.GCA_000347455.1.24.dna_sm.toplevel.fa | sed "s/>/>Plant:Triticum_urartu:/" >> Plant5.fa



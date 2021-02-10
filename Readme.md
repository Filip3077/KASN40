# Detta är en readme 

För att läsa i filer från andra mappar än den som skriptet körs ifrån:

- Om filerna ligger i en undermapp från där skriptet körs: ./<namn på mappen>/filen.msa
- Om filerna ligger i en övermapp till där skriptet körs: ../filen.msa
- Om filerna ligger i en "parralell" mapp till där skriptet körs: ../<namn på mappen>/filen.msa

D.v.s ./ betyder mappen som skriptet är i, ../ betyder "gå upp" en mapp och / är var den ska leta. Så ska man upp flera mappar skriver man ../../ osv ../filen.msa

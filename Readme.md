# En Readme om hur jag tänker att vi kan använda GitHub

Som ni ser har jag "fixat" en del... Jag borde skrivit fixar git för att vara tydligare men det är basicly det jag gjort iaf.

Min tanke är att vi använder git mer som det är skapt för nu genom att man gör sin branch när man vill fixa något, det kan vara nya skript, edits på befintliga, lägga till filer etc. Jag har inte bestämt mig än, men tror det blir bäst att använda mappar som "Filip v.6" tills vi har en konkret skelletmetod som vi alla jobbar i. 

### Arbetsgång:
Eftersom vi alla har GitHub desktop nu (visst?) så är det smidigaste att man gör en branch (har en branch för den veckan?) på hemsidan, då kopieras allt som finns i main till den. Sen går man in i GitHub desktop och "fetchar" sin branch, då finns alla filer och mappar i mappen på ens dator (för mig är det C:\Users\Lucas\Documents\GitHub\KASN40) och man kan köra alla skript och nå alla filer man behöver därifrån. Vi kan lägga in nya spectrum och sånt i Spectra-mappen och göra egna mappar för specifika skript/projekt/test/etc. När man är klar eller vill ladda upp sparar man sina ändringar och går in i GitHub desktop, skriver en liten rubrik/kommentar för vad man gjort och klickar på push origin så läggs det upp i branchen. Sen kan man slå ihop det med main med pull requests. Man kan såklart också jobba direkt i main på sin egen dator, det gjorde jag nu. Jag tror det här kommer funka ganska smidigt och vi kan egentligen arbeta på mer eller mindre samma grejer samtidigt. 

Anledningen till att jag vill ha mappindelningen är att det blir mer städat men då kan det bli så att man får ha massa duplett filer för spectrum som används i flera skript. Som tur är fick jag hälp av min bror hur man hämtar filer från "närliggande" mappar:

#### För att läsa in filer från andra mappar än den som skriptet körs ifrån:

- Om filerna ligger i en undermapp från där skriptet körs: ./<namn på mappen>/filen.msa
- Om filerna ligger i en övermapp till där skriptet körs: ../filen.msa
- Om filerna ligger i en "parralell" mapp till där skriptet körs: ../<namn på mappen>/filen.msa

D.v.s ./ betyder mappen som skriptet är i, ../ betyder "gå upp" en mapp och / är var den ska leta. Så ska man upp flera mappar skriver man ../../ osv. 

T.ex ../../LucasSpectrum/LucPureAg.msa för att gå från KASN40 mappeen i GitHub och hämta en LucPureAg.msa i LucasSpectrum mappen om den ligger i C:\Users\Lucas\Documents ...

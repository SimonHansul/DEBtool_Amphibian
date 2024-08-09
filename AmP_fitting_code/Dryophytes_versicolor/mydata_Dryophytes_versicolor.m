function [data, auxData, metaData, txtData, weights] = mydata_Dryophytes_versicolor
%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Amphibia'; 
metaData.order      = 'Anura'; 
metaData.family     = 'Hylidae';
metaData.species    = 'Dryophytes_versicolor'; 
metaData.species_en = 'Eastern Grey Tree Frog'; 

% unsure about these data, but they are not used for the fit of the parameters
metaData.ecoCode.climate = {'Am', 'BSk', 'BSh'};
metaData.ecoCode.ecozone = {'THn', 'TN'};
metaData.ecoCode.habitat = {'0jFp', 'jiTh', 'jiTg'};
metaData.ecoCode.embryo  = {'Fpf'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biCi'};
metaData.ecoCode.gender  = {'Dg'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(20.5); % K, body temp
metaData.data_0     = {'am'; 'Lp'; 'Li'; 'Wwb'; 'Wwp'; 'Wwi'; 'Ri'}; 
metaData.data_1     = {'L-dL'}; 

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011

metaData.author      = {'Carlo Romoli'};    
metaData.date_subm   = [2024 06 25];              
metaData.email       = {'carlo.romoli@ibacon.com'};            
metaData.address     = {'ibacon GmbH, Germany'};   

%% set data
% zero-variate data

data.am = 7.8*365; units.am = 'd';    label.am = 'life span';                        bibkey.am = 'AnAge';   
  temp.am = C2K(17.4);  units.temp.am = 'K'; label.temp.am = 'temperature';
  comment.am = 'uncertain';

data.ah = 4;  units.ah = 'd';    label.ah = 'age at hatch';                        bibkey.ah = 'AnAge';   
    temp.ah = C2K(17.4);  units.temp.ah = 'K'; label.temp.ah = 'temperature';
    comment.ah = 'uncertain temperature, days to hatching confirmed by Beachy1999';

data.Lp  = 4.56;  units.Lp  = 'cm';  label.Lp  = 'SVL at puberty';                    bibkey.Lp  = 'Höbel2021';
  comment.Lp = '1 year old female';

data.Li  = 6.0;  units.Li  = 'cm';  label.Li  = 'ultimate SVL for females';          bibkey.Li  = 'Wright1949';

% data.Lim  = 5.1; units.Lim  = 'cm'; label.Lim  = 'ultimate SVL for males';           bibkey.Lim  = 'Wright1949';

data.Wwb =  0.0061; units.Wwb = 'g'; label.Wwb = 'wet weight at birth';               bibkey.Wwb = 'Smith2004';
  comment.Wwb = 'at 19 °C';

data.ab = 7;  units.ab = 'd';    label.ab = 'age at birth';                        bibkey.ah = 'AnAge';   
    temp.ab = C2K(27.);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
    comment.ab = 'temperature a bit uncertain';

data.Wwp = 7.2;  units.Wwp = 'g';  label.Wwp = 'wet weight at puberty';              bibkey.Wwp = {'AnAge'};
  comment.Wwp = 'uncertain';
data.Wwi = 11.7; units.Wwi = 'g';   label.Wwi = 'ultimate wet weight for females';    bibkey.Wwi = 'Kuczynski2017';
%data.Wwim = 3.5; units.Wwim = 'g';  label.Wwim = 'ultimate wet weight for males';    bibkey.Wwim = 'HernRami2018';

data.Ri  = 1800/365; units.Ri  = '#/d'; label.Ri  = 'reprod rate';     bibkey.Ri  = 'Höbel2021';   comment.Ri = {"appriximate value considering 2 clutches"};
  temp.Ri  = C2K(20.5);  units.temp.Ri = 'K'; label.temp.Ri = 'temperature';

data.Wg46 = 0.484; units.Wg46 = 'g'; label.Wg46 = 'weight at end of metamorphosis'; bibkey.Wg46 = 'Boone2006'; comment.Wg46 = 'embryos at 24°C till G25, then at outside temp';
data.ag46 = 25.4; units.ag46 = 'd'; label.ag46 = 'time to end of metamorphosis from birth'; bibkey.ag46 = 'Boone2006'; comment.ag46 = 'embryos at 24°C till G25, then at outside temp';

data.W2g46 = 0.56; units.W2g46 = 'g'; label.W2g46 = 'weight at end of metamorphosis'; bibkey.W2g46 = 'Macklem2013'; comment.W2g46 = 'at 23°C';
data.a2g46 = 35.5; units.a2g46 = 'd'; label.a2g46 = 'time to end of metamorphosis from birth'; bibkey.a2g46 = 'Macklem2013'; comment.a2g46 = 'at 23°C';



%% uni-variate data

% Smith 2004 (initially approx. Average value of low density with intraspecific and interspecific competition)
data.Ww6w = [ ... % time (d), mass (g)
     0    0.0061
     42   0.09
];
units.Ww6w = {'d', 'g'}; label.Ww6w = {'time since birth', 'wet mass (g)'};
temp.Ww6w = C2K(19); units.temp.Ww6w = 'K'; label.temp.Ww6w = 'temperature';
bibkey.Ww6w = 'Smith2004'; comment.Ww6w = 'data for tadpoles';

data.LLL = [ ... % time (d) , mass (g)
10.0  0.041176471
20.0  0.179705882
30.0  0.353529412
40.0  0.422352941
50.72	0.310294118
];
units.LLL = {'d', 'g'}; label.LLL = {'time since hatch', 'wet mass (g)'};
temp.LLL = C2K(22); units.temp.LLL = 'K'; label.temp.LLL = 'temperature';
bibkey.LLL = 'Beachy1999'; comment.LLL = 'data for tadpoles';

data.LHH = [ ... % time (d) , mass (g)
10.0	0.040294118
20.0	0.209705882
30.0	0.504411765
40.0	0.613823529
42.78	0.430294118
];
units.LHH = {'d', 'g'}; label.LHH = {'time since hatch', 'wet mass (g)'};
temp.LHH = C2K(22); units.temp.LHH = 'K'; label.temp.LHH = 'temperature';
bibkey.LHH = 'Beachy1999'; comment.LHH = 'data for tadpoles';

data.LLH = [ ... % time (d) , mass (g)
10.0	0.040294118
20.0	0.178823529
30.0	0.347352941
40.0	0.527352941
46.89	0.433823529
];
units.LLH = {'d', 'g'}; label.LLH = {'time since hatch', 'wet mass (g)'};
temp.LLH = C2K(22); units.temp.LLH = 'K'; label.temp.LLH = 'temperature';
bibkey.LLH = 'Beachy1999'; comment.LLH = 'data for tadpoles';

data.HHH = [ ... % time (d) , mass (g)
10.0	0.03802303
20.0	0.226634425
30.0	0.52707348
40.0	0.566569493
42.42	0.481470933
];
units.HHH = {'d', 'g'}; label.HHH = {'time since hatch', 'wet mass (g)'};
temp.HHH = C2K(22); units.temp.HHH = 'K'; label.temp.HHH = 'temperature';
bibkey.HHH = 'Beachy1999'; comment.HHH = 'data for tadpoles';

data.HLL = [ ... % time (d) , mass (g)
10.0	0.038020663
20.0	0.224856903
30.0	0.401037479
40.0	0.45562231
47.84	0.306859593
];
units.HLL = {'d', 'g'}; label.HLL = {'time since hatch', 'wet mass (g)'};
temp.HLL = C2K(22); units.temp.HLL = 'K'; label.temp.HLL = 'temperature';
bibkey.HLL = 'Beachy1999'; comment.HLL = 'data for tadpoles';

data.HHL = [ ... % time (d) , mass (g)
10.0	0.03890824
20.0	0.22574448
30.0	0.477369141
40.0	0.45562231
41.98	0.371392392
];
units.HHL = {'d', 'g'}; label.HHL = {'time since hatch', 'wet mass (g)'};
temp.HHL = C2K(22); units.temp.HHL = 'K'; label.temp.HHL = 'temperature';
bibkey.HHL = 'Beachy1999'; comment.HHL = 'data for tadpoles';



%% data from our experiment
data.temptimeH = [ ... % temp(C), time(d)
    18 41.7
    20 36.7
    22 34.8
    24 30
];
units.temptimeH = {'°C', 'g'}; label.temptimeH = {'temperature', 'time (d)'};
bibkey.temptimeH = 'Experiment'; comment.temptimeH = 'data for tadpoles';

data.temptimeH2 = [ ... % temp(C), time(d)
    22 40
    24 36
];
units.temptimeH2 = {'°C', 'g'}; label.temptimeH2 = {'temperature', 'time (d)'};
bibkey.temptimeH2 = 'Experiment w starv 6w'; comment.temptimeH2 = 'data for tadpoles';

data.temptimeL = [ ... % temp(C), time(d)
    18 62.2
    20 49.2
    22 41.7
    24 40
];
units.temptimeL = {'°C', 'g'}; label.temptimeL = {'temperature', 'time (d)'};
bibkey.temptimeL = 'Experiment'; comment.temptimeL = 'data for tadpoles';

data.tempDwH = [... % temp(C), dry weight(g)
    18 0.0716
    20 0.0742
    22 0.0517
    24 0.0546
];
units.tempDwH = {'°C', 'g'}; label.tempDwH = {'temperature', 'dry weight (g)'};
bibkey.tempDwH = 'Experiment'; comment.tempDwH = 'data for tadpoles';

data.tempDwH2 = [... % temp(C), dry weight(g)
    22 0.0408
    24 0.0417
];
units.tempDwH2 = {'°C', 'g'}; label.tempDwH2 = {'temperature', 'dry weight (g)'};
bibkey.tempDwH2 = 'Experiment w starv 6w'; comment.tempDwH2 = 'data for tadpoles';

data.tempDwL = [... % temp(C), dry weight(g)
    18 0.0228
    20 0.0216
    22 0.0236
    24 0.0175
];
units.tempDwL = {'°C', 'g'}; label.tempDwL = {'temperature', 'dry weight (g)'};
bibkey.tempDwL = 'Experiment'; comment.tempDwL = 'data for tadpoles';

data.timeWwH18 = [... % time(d), wet weight(g)
    7	0.0915
    13	0.06495
    19	0.24445
    26	0.10165
    33	0.2084
];
units.timeWwH18 = {'d', 'g'}; label.timeWwH18 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwH18 = 'Experiment weekly measurments'; comment.timeWwH18 = 'data for tadpoles';
temp.timeWwH18 = C2K(18);  units.temp.timeWwH18 = 'K'; label.temp.timeWwH18 = 'temperature';

data.timeWwH20 = [... % time(d), wet weight(g)
    7	0.046
    13	0.0654
    20	0.09035
    29	0.1069
    34	0.47625
];
units.timeWwH20 = {'d', 'g'}; label.timeWwH20 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwH20 = 'Experiment weekly measurments'; comment.timeWwH20 = 'data for tadpoles';
temp.timeWwH20 = C2K(20);  units.temp.timeWwH20 = 'K'; label.temp.timeWwH20 = 'temperature';

data.timeWwH22 = [... % time(d), wet weight(g)
    7	0.056166667
    14	0.155533333
    21	0.1539
    30	0.261733333
    37  0.3679
];
units.timeWwH22 = {'d', 'g'}; label.timeWwH22 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwH22 = 'Experiment weekly measurments'; comment.timeWwH22 = 'data for tadpoles';
temp.timeWwH22 = C2K(22);  units.temp.timeWwH22 = 'K'; label.temp.timeWwH22 = 'temperature';

data.timeWwH24 = [... % time(d), wet weight(g)
    7	0.0685
    15	0.0612
    22	0.0576
];
units.timeWwH24 = {'d', 'g'}; label.timeWwH24 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwH24 = 'Experiment weekly measurments'; comment.timeWwH24 = 'data for tadpoles';
temp.timeWwH24 = C2K(24);  units.temp.timeWwH24 = 'K'; label.temp.timeWwH24 = 'temperature';




% low food level

data.timeWwL18 = [... % time(d), wet weight(g)
7	0.089
13	0.08005
19	0.08585
26	0.09235
33	0.1114
44	0.25185
51	0.8075
58	0.339
];
units.timeWwL18 = {'d', 'g'}; label.timeWwL18 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwL18 = 'Experiment weekly measurments'; comment.timeWwL18 = 'data for tadpoles';
temp.timeWwL18 = C2K(18);  units.temp.timeWwL18 = 'K'; label.temp.timeWwL18 = 'temperature';

data.timeWwL20 = [... % time(d), wet weight(g)
7	0.02825
13	0.09755
20	0.12465
29	0.08055
34	0.12175
44	0.4186
];
units.timeWwL20 = {'d', 'g'}; label.timeWwL20 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwL20 = 'Experiment weekly measurments'; comment.timeWwL20 = 'data for tadpoles';
temp.timeWwL20 = C2K(20);  units.temp.timeWwL20 = 'K'; label.temp.timeWwL20 = 'temperature';

data.timeWwL22 = [... % time(d), wet weight(g)
7	0.09675
14	0.13035
21	0.0948
30	0.19595
37	0.40155
];
units.timeWwL22 = {'d', 'g'}; label.timeWwL22 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwL22 = 'Experiment weekly measurments'; comment.timeWwL22 = 'data for tadpoles';
temp.timeWwL22 = C2K(22);  units.temp.timeWwL22 = 'K'; label.temp.timeWwL22 = 'temperature';

data.timeWwL24 = [... % time(d), wet weight(g)
7	0.04925
15	0.0925
22	0.1043
30	0.15855
37	0.20215
];
units.timeWwL24 = {'d', 'g'}; label.timeWwL24 = {'time (d)', 'wet weight (g)'};
bibkey.timeWwL24 = 'Experiment weekly measurments'; comment.timeWwL24 = 'data for tadpoles';
temp.timeWwL24 = C2K(24);  units.temp.timeWwL24 = 'K'; label.temp.timeWwL24 = 'temperature';



%% set weights for all real data
weights = setweights(data, []);
weights.Li = 5 * weights.Li;
weights.temptimeL = repelem(1,length(data.temptimeL))';
weights.temptimeH = repelem(1,length(data.temptimeH))';
weights.temptimeH2 = repelem(1,length(data.temptimeH2))';

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Group plots
set1 = {'timeWwH18','timeWwH20', 'timeWwH22', 'timeWwH24'}; subtitle1 = {'Different temperature at high food'};
set2 = {'timeWwL18', 'timeWwL20','timeWwL22', 'timeWwL24'}; subtitle2 = {'Different temperature at low food'};
set3 = {'LLL','LHH', 'LLH', 'HHH', 'HLL', 'HHL'}; subtitle3 = {'Different feeding regimes'};
metaData.grp.sets = {set1, set2,set3};
metaData.grp.subtitle = {subtitle1, subtitle2, subtitle3};

%% Discussion points
% D1 = 'Males are assumed to differ from females by {p_Am} only';
% D2 = 'Temperatures are guessed';
% metaData.discussion = struct('D1', D1, 'D2', D2);

%% Links
metaData.links.id_CoL = '37VNH'; % Cat of Life
metaData.links.id_ITIS = '1104993'; % ITIS
metaData.links.id_EoL = '1019447'; % Ency of Life
metaData.links.id_Wiki = 'Dryophytes_eximius'; % Wikipedia
metaData.links.id_ADW = 'Hyla_eximia'; % ADW
metaData.links.id_Taxo = '138153'; % Taxonomicon
metaData.links.id_WoRMS = ''; % WoRMS
metaData.links.id_amphweb = 'Hyla+eximia'; % AmphibiaWeb


%% References
bibkey = 'EoL'; type = 'Misc'; bib = ...
'howpublished = {\url{https://eol.org/pages/52234580}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'AnAge'; type = 'Misc'; bib = ...
'howpublished = {\url{https://genomics.senescence.info/species/entry.php?species=Hyla_versicolor}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Hoebel2021'; type = 'Article'; bib = [ ...  
'doi = {10.1163/15685381-bja10077}, ' ...
'author = {Gerlinde Höbel and Robb Kolodziej and Dustin Nelson and Christopher White}, ' ... 
'year = {2021}, ' ...
'title = {Effect of body size, age and timing of breeding on clutch and egg size of female Eastern Gray Treefrogs, Hyla versicolor}, ' ...
'journal = {Amphibia-Reptilia}, ' ...
'volume = {43}, ' ...
'number = {1},' ...
'pages = {25-35}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Wright1949'; type = 'Book'; bib = [ ... 
'author = {Wright, Albert Hazen and Wright, Anna Allen}, ' ...
'year = {1949}, ' ...
'title  = {Handbook of Frogs and Toads of the United States and Canada}, ' ...
'publisher = {Cornell Univ. Press}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Smith04'; type = 'Article'; bib = [ ...  
'doi = {https://doi.org/10.1111/j.0030-1299.2004.12841.x}, ' ...
'author = {Smith, Geoffrey R. and Dingfelder, Haley A. and Vaala, David A.}, ' ... 
'year = {2004}, ' ...
'title = {Asymmetric competition between Rana clamitans and Hyla versicolor tadpoles},' ...
'journal = {Oikos},' ...
'volume = {105}, ' ...
'number = {3},' ...
'pages = {626-632}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kuczynski2017'; type = 'Article'; bib = [ ...  
'doi = {https://doi.org/10.1016/j.beproc.2016.11.019},' ...
'author = {Michael C. Kuczynski and Thomas Getty and Eben Gering},' ... 
'year = {2017}, ' ...
'title = {Larger females are choosier in the gray treefrog (Hyla versicolor)},' ...
'journal = {Behavioural Processes},' ...
'volume = {135}, ' ...
'pages = {29-35}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Boone2006'; type = 'Article'; bib = [ ...  
'doi = {https://doi.org/10.1897/06-235R.1},' ...
'author = {Boone, Michelle D. and Bridges-Britton, Christine M.},' ... 
'year = {2006}, ' ...
'title = {Examining multiple sublethal contaminants on the gray treefrog (Hyla versicolor): Effects of an insecticide, herbicide, and fertilizer},' ...
'journal = {Environmental Toxicology and Chemistry},' ...
'volume = {25}, ' ...
'number = {12},' ...
'pages = {3261-3265}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Macklem2013'; type = 'misc'; bib = [ ...  
'url = {https://digitalcommons.lib.uconn.edu/cgi/viewcontent.cgi?article=1012&context=srhonors_holster},' ...
'author = {Diana Cristina Macklem},' ... 
'year = {2013}, ' ...
'title = {The Effects of Temperature Frequency and Magnitude on Lithobates sylvaticus and Hyla versicolor},'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Beachy99'; type = 'Article'; bib = [ ...  
'doi = {https://doi.org/10.1002/(SICI)1097-010X(19990501)283:6<522::AID-JEZ3>3.0.CO;2-3},' ...
'author = {Beachy, Christopher K. and Surges, Tammy H. and Reyes, Monica},' ... 
'year = {1999}, ' ...
'title = {Effects of developmental and growth history on metamorphosis in the gray treefrog, Hyla versicolor (Amphibia, Anura)},' ...
'journal = {Journal of Experimental Zoology},' ...
'volume = {283}, ' ...
'number = {6},' ...
'pages = {522-530}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
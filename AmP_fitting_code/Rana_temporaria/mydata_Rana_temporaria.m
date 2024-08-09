function [data, auxData, metaData, txtData, weights] = mydata_Rana_temporaria

%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Amphibia'; 
metaData.order      = 'Anura'; 
metaData.family     = 'Ranidae';
metaData.species    = 'Rana_temporaria'; 
metaData.species_en = 'Common frog'; 

metaData.ecoCode.climate = {'Cfb', 'Dfb', 'Dfc'};
metaData.ecoCode.ecozone = {'THp'};
metaData.ecoCode.habitat = {'0jFp', 'jiTf', 'jiTg', 'jiTi'};
metaData.ecoCode.embryo  = {'Fpf'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'bjO', 'jiCi'};
metaData.ecoCode.gender  = {'Dg'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(18); % K, body temp
metaData.data_0     = {'ab'; 'aj'; 'ap'; 'am'; 'Lj'; 'Li'; 'Wwb'; 'Wwi'; 'Ri'}; 
metaData.data_1     = {'L-dL'}; 

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011

metaData.author   = {'Bas Kooijman'};    
metaData.date_subm = [2017 10 28];              
metaData.email    = {'bas.kooijman@vu.nl'};            
metaData.address  = {'VU University, Amsterdam'};   

metaData.curator     = {'Starrlight Augustine'};
metaData.email_cur   = {'starrlight@akvaplan.niva.no'}; 
metaData.date_acc    = [2017 10 28]; 


%% set data
% zero-variate data

data.ab = 35;      units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = 'ADW';   
  temp.ab = C2K(10);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
  comment.ab = 'breeding in early spring';
data.tj = 1.5*30.5;  units.tj = 'd';    label.tj = 'time since birth at metam'; bibkey.tj = 'amphibiaweb';
  temp.tj = C2K(18);  units.temp.tj = 'K'; label.temp.tj = 'temperature';
data.tp = 1095;  units.tp = 'd';    label.tp = 'time since metam at puberty'; bibkey.tp = 'AnAge';
  temp.tp = C2K(18);  units.temp.tp = 'K'; label.temp.tp = 'temperature';
data.am = 14*365;  units.am = 'd';    label.am = 'life span';                bibkey.am = 'AnAge';   
  temp.am = C2K(18);  units.temp.am = 'K'; label.temp.am = 'temperature'; 

data.Lj  = 1.5;     units.Lj  = 'cm';  label.Lj  = 'body length at metam'; bibkey.Lj  = 'kikkersite'; 
data.Li  = 9;    units.Li  = 'cm';  label.Li  = 'ultimate body length';   bibkey.Li  = 'Wiki';

data.Wwb = 0.1;  units.Wwb = 'g';   label.Wwb = 'wet weight at birth';     bibkey.Wwb = 'Ruthsatz20a';
  comment.Wwb = 'from paper';
data.Wwi = 68.6;  units.Wwi = 'g';   label.Wwi = 'ultimate wet weight';     bibkey.Wwi = 'Wiki';
  comment.Wwi = 'weight at 9 cm for Rana lessonae: (9/7.6)^3*41.3';

data.Ri  = 1500/365; units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';  bibkey.Ri  = 'Wiki';   
  temp.Ri = C2K(18);  units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
  comment.Ri = '1000-2000 eggs per season';

  
% % univariate data
% % length-change in length
data.LdL = [... % body length (cm), change in length (cm/76 d)
3.135	1.760
3.485	1.913
3.546	1.485
3.583	1.934
3.649	1.782
3.658	2.148
3.705	1.843
3.798	1.406
3.821	1.175
3.859	1.559
3.877	1.424
3.891	1.537
3.896	1.188
3.901	1.934
3.915	2.052
3.915	1.755
3.924	1.445
3.961	2.039
3.971	1.275
3.980	1.952
3.980	1.428
4.004	1.930
4.008	1.384
4.027	1.590
4.036	1.760
4.060	1.808
4.064	1.620
4.102	1.376
4.139	1.873
4.153	1.624
4.158	1.716
4.158	1.371
4.158	1.358
4.172	1.092
4.176	1.039
4.204	1.882
4.218	1.725
4.218	1.445
4.223	1.786
4.251	1.550
4.251	1.581
4.256	1.655
4.260	1.463
4.260	1.271
4.279	1.100
4.293	1.624
4.298	1.716
4.316	1.817
4.316	1.306
4.321	1.550
4.340	1.131
4.344	1.900
4.344	1.410
4.368	1.559
4.382	1.725
4.410	1.441
4.438	1.507
4.447	1.773
4.456	1.668
4.466	1.245
4.466	1.197
4.480	1.515
4.489	1.734
4.494	1.410
4.498	1.148
4.498	0.991
4.503	1.114
4.526	1.066
4.526	1.009
4.550	1.563
4.554	1.197
4.564	0.956
4.573	1.122
4.582	1.742
4.587	0.904
4.601	1.467
4.624	1.459
4.624	0.808
4.662	1.528
4.662	1.410
4.662	1.131
4.704	0.991
4.718	0.904
4.736	1.031
4.741	1.162
4.755	1.096
4.774	0.847
4.792	1.131
4.816	1.520
4.816	1.515
4.825	0.974
4.834	1.109
4.844	1.328
4.848	1.288
4.909	1.306
4.909	0.769
4.970	1.140
4.988	0.817
4.993	1.205
5.026	0.729
5.030	0.943
5.044	0.659
5.152	0.917
5.161	1.467
5.170	1.441
5.203	0.991
5.203	1.070
5.217	0.651
5.287	0.817
5.651	0.747
6.197	0.358];
data.LdL(:,2) = data.LdL(:,2)/ 76; % convert cm/76 d to cm/d
units.LdL  = {'cm', 'cm/d'}; label.LdL = {'body length', 'change in body length'};  
temp.LdL   = C2K(18); units.temp.LdL = 'K'; label.temp.LdL = 'temperature';
bibkey.LdL = 'LomaLard2009';



%% additional data Ruthsatz
data.SVLd = [... % age since hatch (d), SVL (mm)
13.66	0.677
16.00	0.756
18.05	0.914
19.68	1.061
20.96	1.05
28.58	1.082
35.82	1.141
];
units.SVLd  = {'d', 'cm'}; label.SVLd = {'age from birth', 'SVL'};  
temp.SVLd   = C2K(22); units.temp.SVLd = 'K'; label.temp.SVLd = 'temperature';
bibkey.SVLd = 'Ruthsatz20a';



data.massd = [... % age since hatch (d), wet weight (mm)
13.66	0.09702
16.00	0.13673
18.05	0.28161
19.68	0.40154
20.96	0.362
28.58	0.10521 
35.82	0.14313
];
units.massd  = {'d', 'g'}; label.massd = {'age from birth', 'wet weight'};  
temp.massd   = C2K(22); units.temp.massd = 'K'; label.temp.massd = 'temperature';
bibkey.massd = 'Ruthsatz20a';

data.dstage = [... % age since hatch (d), SVL (mm)
25  13.66
42  20.96
46  28.58
];
units.dstage  = {'G', 'd'}; label.dstage = {'stage', 'time since hatch'};  
temp.dstage   = C2K(22); units.temp.dstage = 'K'; label.temp.dstage = 'temperature';
bibkey.dstage = 'Ruthsatz20a';


data.tbirth=13.66; units.tbirth = 'd';    label.tbirth = 'age at G25';   bibkey.tbirth = 'Ruthsatz20a';   
  temp.tbirth = C2K(22);  units.temp.tbirth = 'K'; label.temp.tbirth = 'temperature';
data.tmetsta=20.96; units.tmetsta = 'd';    label.tmetsta = 'age at G42';   bibkey.tmetsta = 'Ruthsatz20a';   
  temp.tmetsta = C2K(22);  units.temp.tmetsta = 'K'; label.temp.tmetsta = 'temperature';
data.tmetsto=28.58; units.tmetsto = 'd';    label.tmetsto = 'age at G46';   bibkey.tmetsto = 'Ruthsatz20a';   
  temp.tmetsto = C2K(22);  units.temp.tmetsto = 'K'; label.temp.tmetsto = 'temperature';


% relative values
data.dtemp = [... % age since hatch (d), SVL (mm)
14  35.8
18  28.5
22  20.9
25  15  
28  13.1
];
units.dtemp  = {'°C', 'd'}; label.dtemp = {'temperature', 'duration G42'};  
bibkey.dtemp = 'Ruthsatz20b';


data.dtotltemp = [... % age since hatch (d), total length (mm)
14  4.3
18  4.1
22  3.0
25  2.2
28  1.9
];
units.dtotltemp  = {'°C', 'cm'}; label.dtotltemp = {'temperature', 'total length G42'};  
bibkey.dtotltemp = 'Ruthsatz20b';

data.dsvlltemp = [... % age since hatch (d), total length (mm)
14  1.6
18  1.6
22  1.1
25  0.9
28  0.8
];
units.dsvlltemp  = {'°C', 'cm'}; label.dsvlltemp = {'temperature', 'SVL length G42'};  
bibkey.dsvlltemp = 'Ruthsatz20b';

data.dmasstemp = [... % age since hatch (d), mass (g)
14  0.603
18  0.483
22  0.361
25  0.131
28  0.083
];
units.dmasstemp  = {'°C', 'g'}; label.dmasstemp = {'temperature', 'wet weight G42'};  
bibkey.dmasstemp = 'Ruthsatz20b';

%% data from Morand97
data.tempdry = [... % temperature, dry weight at end of metamorphosis (mg)
    15    0.041
    27    0.024
];
units.tempdry  = {'°C', 'g'}; label.tempdry = {'temperature', 'dry weight G46'};  
bibkey.tempdry = 'Morand97';

data.tempsdur = [... % temperature, dry weight at end of metamorphosis (mg)
    15    64.26
    27    48.29
];
units.tempsdur  = {'°C', 'd'}; label.tempsdur = {'temperature', 'days G46'};  
bibkey.tempsdur = 'Morand97';

%% data from Aviles22
% additional point to pin the metamorphosis timing
data.startj=48; units.startj = 'd';    label.startj = 'age at G42';  bibkey.startj = 'Aviles22';   
  temp.startj = C2K(16.5);  units.temp.startj = 'K'; label.temp.startj = 'temperature';

data.ww42=0.7855; units.ww42 = 'g'; label.ww42 = 'wet weight at G42'; bibkey.ww42 = 'Aviles22';
  temp.ww42 = C2K(16.5); units.temp.ww42 = 'K'; label.temp.ww42 = 'temperature';

data.svl42=1.8; units.svl42 = 'cm'; label.svl42 = 'SVL at G42'; bibkey.svl42 = 'Aviles22';
  temp.svl42 = C2K(16.5); units.temp.svl42 = 'K'; label.temp.svl42 = 'temperature';

data.totlen42=4.4; units.totlen42 = 'cm'; label.totlen42 = 'total lenght at G42'; bibkey.totlen42 = 'Aviles22';
  temp.totlen42 = C2K(16.5); units.temp.totlen42 = 'K'; label.temp.totlen42 = 'temperature';


%% set weights for all real data
weights = setweights(data, []);

weights.dstage=[5;5;5];  % artificially increase the weight of the timing information so that the model is forced to pass thorough the developmental stages at the right moment.
%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Facts
F1 = 'Males slightly smaller than females';
metaData.bibkey.F1 = 'Wiki'; 
metaData.facts = struct('F1',F1);

%% Links
metaData.links.id_CoL = '6WR99'; % Cat of Life
metaData.links.id_ITIS = '173444'; % ITIS
metaData.links.id_EoL = '331133'; % Ency of Life
metaData.links.id_Wiki = 'Rana_temporaria'; % Wikipedia
metaData.links.id_ADW = 'Rana_temporaria'; % ADW
metaData.links.id_Taxo = '47660'; % Taxonomicon
metaData.links.id_WoRMS = '1488291'; % WoRMS
metaData.links.id_amphweb = 'Rana+temporaria'; % AmphibiaWeb
metaData.links.id_AnAge = 'Rana_temporaria'; % AnAge



%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Rana_temporaria}}';
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
bibkey = 'LomaLard2009'; type = 'Article'; bib = [ ... 
'author = {J. Loman and B. Lardner}, ' ... 
'year = {2009}, ' ...
'title = {Density dependent growth in adult brown frogs \emph{Rana arvalis} and \emph{Rana temporaria} - A field experiment}, ' ...
'journal = {Acta Oecologica}, ' ...
'volume = {35}, ' ...
'pages = {824-830}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'amphibiaweb'; type = 'Misc'; bib = ...
'howpublished = {\url{http://amphibiaweb.org/cgi/amphib_query?rel-common_name=like&where-scientific_name=Rana+temporaria}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'AnAge'; type = 'Misc'; bib = ...
'howpublished = {\url{http://genomics.senescence.info/species/entry.php?species=Rana_temporaria}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'ADW'; type = 'Misc'; bib = ...
'howpublished = {\url{http://animaldiversity.org/accounts/Rana_temporaria/}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'kikkersite'; type = 'Misc'; bib = ...
'howpublished = {\url{https://www.kikkersite.nl/bruine.php}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Aviles22'; type = 'Article'; bib = [ ... 
'author = {Amandine Aviles and Katrine Hulgard and John W. Green and Annette Duus and Bente Holbech and Jane E. Morthorst},' ... 
'year = {2022}, ' ...
'title = {Effects of sodium perchlorate and 6-propylthiouracil on metamorphosis and thyroid gland histopathology in the European common frog (Rana temporaria)},' ...
'journal = {Aquatic Toxicology}, ' ...
'volume = {244}, ' ...
'pages = 106094},' ...
'doi = {https://doi.org/10.1016/j.aquatox.2022.106094}'
];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
%
bibkey = 'Ruthsatz20b'; type = 'Article'; bib = [ ... 
'author = {Ruthsatz, Katharina and Dausmann, Kathrin H and Paesler, Katharina and Babos, Patricia and Sabatino, Nikita M and Peck, Myron A and Glos, Julian}, ' ... 
'year = {2020}, ' ...
'title = {Shifts in sensitivity of amphibian metamorphosis to endocrine disruption: the common frog (Rana temporaria) as a case study},' ...
'journal = {Conservation Physiology}, ' ...
'volume = {8},' ...
'number = {1},' ...
'pages = {coaa100},' ...
'doi = {10.1093/conphys/coaa100}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
%
bibkey = 'Ruthsatz20a'; type = 'Article'; bib = [ ... 
'author = {Katharina Ruthsatz and Kathrin H. Dausmann and Claudia Drees and Laura I. Becker and Lisa Hartmann and Janica Reese and Steffen Reinhardt and Tom Robinson and Nikita M. Sabatino and Myron A. Peck and Julian Glos}, ' ... 
'year = {2020}, ' ...
'title = {Altered thyroid hormone levels affect the capacity for temperature-induced developmental plasticity in larvae of Rana temporaria and Xenopus laevis}, ' ...
'journal = {Journal of Thermal Biology}, ' ...
'volume = {90}, ' ...
'pages = {102599},'...
'doi = {https://doi.org/10.1016/j.jtherbio.2020.102599},'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
%
bibkey = 'LomaLard2009'; type = 'Article'; bib = [ ... 
'author = {J. Loman and B. Lardner}, ' ... 
'year = {2009}, ' ...
'title = {Density dependent growth in adult brown frogs \emph{Rana arvalis} and \emph{Rana temporaria} - A field experiment}, ' ...
'journal = {Acta Oecologica}, ' ...
'volume = {35}, ' ...
'pages = {824-830}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Morand97'; type = 'Article'; bib = [ ... 
'author = {Alain Morand and Pierre Joly and Odile Grolet}, ' ... 
'year = {1997}, ' ...
'title = {Phenotypic variation in metamorphosis in five anuran species along a gradient of stream influence}, ' ...
'journal = {Comptes Rendus De L Academie Des Sciences Serie Iii-sciences De La Vie-life Sciences}, ' ...
'volume={320}, ' ...
'pages={645-652}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];



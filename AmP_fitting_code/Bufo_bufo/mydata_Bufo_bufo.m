function [data, auxData, metaData, txtData, weights] = mydata_Bufo_bufo

%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Amphibia'; 
metaData.order      = 'Anura'; 
metaData.family     = 'Bufonidae';
metaData.species    = 'Bufo_bufo'; 
metaData.species_en = 'Common toad'; 

metaData.ecoCode.climate = {'Cfb', 'Dfb', 'Dfc'};
metaData.ecoCode.ecozone = {'THp'};
metaData.ecoCode.habitat = {'0jFp', 'jiTf', 'jiTg'};
metaData.ecoCode.embryo  = {'Fh'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'bjO', 'jiCi'};
metaData.ecoCode.gender  = {'Dg'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(15); % K, body temp
metaData.data_0     = {'ab'; 'aj'; 'ap'; 'am'; 'Lj'; 'Lp'; 'Li'; 'Wwb'; 'Wwj'; 'Wwi'; 'Ri'}; 
metaData.data_1     = {'t-L'}; 

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011

metaData.author   = {'Jan Baas';'Bas Kooijman'};    
metaData.date_subm = [2014 05 31];              
metaData.email    = {'janbaa@ceh.ac.uk'};            
metaData.address  = {'VU University, Amsterdam'};   

metaData.author_mod_1   = {'Bas Kooijman'};    
metaData.date_mod_1 = [2017 10 24];              
metaData.email_mod_1    = {'bas.kooijman@vu.nl'};            
metaData.address_mod_1  = {'VU University, Amsterdam'};   

metaData.author_mod_2   = {'Bas Kooijman'};    
metaData.date_mod_2 = [2019 05 15];              
metaData.email_mod_2    = {'bas.kooijman@vu.nl'};            
metaData.address_mod_2  = {'VU University, Amsterdam'};   

metaData.curator     = {'Starrlight Augustine'};
metaData.email_cur   = {'starrlight@akvaplan.niva.no'}; 
metaData.date_acc    = [2019 05 15]; 


%% set data
% zero-variate data

data.ab = 25;    units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = 'KatzWari2003';   
  temp.ab = C2K(15);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
  comment.ab = 'average temperature between 12 and 18 degrees';
data.tj = 58;    units.tj = 'd';    label.tj = 'time since birth at metam';bibkey.tj = 'KatzWari2003';   
  temp.tj = C2K(20.5);  units.temp.tj = 'K'; label.temp.tj = 'temperature';
data.tp = 2.5*365;  units.tp = 'd';    label.tp = 'time since metam at puberty'; bibkey.tp = 'EoL';
  temp.tp = C2K(15);  units.temp.tp = 'K'; label.temp.tp = 'temperature';
data.am = 13000; units.am = 'd';    label.am = 'life span';                bibkey.am = 'EoL';   
  temp.am = C2K(15);  units.temp.am = 'K'; label.temp.am = 'temperature'; comment.am="Maximum lifespan in captivity";

data.Lj  = 1.3;  units.Lj  = 'cm';  label.Lj  = 'total length at metam';   bibkey.Lj  = 'KatzWari2003'; 
data.Lp  = 8.57;   units.Lp  = 'cm';  label.Lp  = 'total length at puberty'; bibkey.Lp  = 'Reading91'; 
data.Li  = 10.8;   units.Li  = 'cm';  label.Li  = 'ultimate total length for females';   bibkey.Li  = 'Hemelaar88';
    comment.Li="France pop. value";
data.Lim  = 8.5;   units.Lim  = 'cm';  label.Lim  = 'ultimate total length for males';   bibkey.Lim  = 'Hemelaar88';
    comment.Lim="France pop. value";

% data.Wwb = 0.0042; units.Wwb = 'g';   label.Wwb = 'wet weight at birth';     bibkey.Wwb = 'amphibiaweb';
%   comment.Wwb = 'based on Anaxyrus americanus';
data.Wwj = 0.2; units.Wwj = 'g';   label.Wwj = 'wet weight at end metam';     bibkey.Wwj = 'KatzWari2003';
data.Wwi = 145;   units.Wwi = 'g';   label.Wwi = 'ultimate wet weight';     bibkey.Wwi = {'KatzWari2003','ReadClar1995'};
  comment.Wwi = 'based on Wwi = exp(-10.1 + 3.22 * 1og(10*Li)), see F4';
data.Wwim = 53.8;   units.Wwim = 'g';   label.Wwim = 'ultimate wet weight of males'; bibkey.Wwim = 'ReadClar1995';
  comment.Wwim = 'based on Wwim = exp(-8.01 + 2.70 * 1og(10*Lim))';
  
%data.Ri  = 1500/365;   units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = 'EoL';   
data.Ri  = 7000/365;   units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = 'Cvetkovic08';   
temp.Ri = C2K(15);  units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
comment.Ri = '1 brood/year; refers to serbian pop. which should have similar size to France pop. used in Hemelaar88';

%% Katzmann data
% after assessing the original paper, the average temperature of the
% experiment was 20.5 degrees
% uni-variate data
% time-length
data.tL = [ ... % time since hatch (day 1 = hatching day), body length (mm)
16.795	3.036
16.796	3.993
23.102	4.949
23.931	6.065
23.932	6.862
30.897	8.855
30.898	9.891
30.899	11.087
37.694	11.884
38.033	9.891
38.034	10.848
44.319	10.928
44.320	11.725
44.321	12.841
50.945	9.971
50.946	11.007
50.947	11.964
50.948	13.000
58.082	11.007
58.083	12.043
58.084	13.000
];
data.tL(:,1) = data.tL(:,1) - 16; % set birth at zero, so easier to analyse
data.tL(:,2) = data.tL(:,2)/10; % convert mm to cm
units.tL   = {'d', 'cm'};  label.tL = {'time since birth', 'SVL'};  
temp.tL =20.5; units.temp.tL = 'C'; label.temp.tL = 'temperature';
%temp.tL    = [0 10 20 30 40 50; 5 10 15 26 26 10]';  units.temp.tL = {'t','C'}; label.temp.tL = {'time','temperature'};
bibkey.tL = 'KatzWari2003';
%comment.tL = 'temperatures ranged from 5 to 26 C; profile is guessed';

data.SL = [ ... % stage (G), SVL (mm)
    25    3.4
    42    11.4
    46    11.9
    ];
units.SL   = {'G', 'mm'};  label.SL = {'stage', 'SVL'};
temp.SL = 20.5; units.temp.SL = 'C'; label.temp.SL = 'temperature';
bibkey.SL = 'KatzWari2003';

data.SW = [ ... % stage (G), Ww (g)
    42    0.35
    46    0.2
    ];
units.SW   = {'G', 'g'};  label.SW = {'stage', 'Wet weight'};
temp.SW = 20.5; units.temp.SW = 'C'; label.temp.SW = 'temperature';
bibkey.SW = 'KatzWari2003';

%% Brunelli data
% experiments start from birth
data.tW = [ ...
    1	0.060759562
    9	0.104913515
    17	0.197144123
    25	0.294374731
    33	0.300836377
    41	0.1992211
    ];
data.tW(:,1) = data.tW(:,1)-1; % the 0 is birth
units.tW   = {'d', 'g'};  label.tW = {'time', 'Wet weight'};
temp.tW = 22; units.temp.tW = 'C'; label.temp.tW = 'temperature';
bibkey.tW = 'BrunBern2009';

data.St = [...
    42 37
    46 48
    ];
units.St = {'G', 'd'};  label.St = {'stage', 'time since birth'};
bibkey.St = 'BrunBern2009';

%% Morand data
% data from Morand97
data.tempdry = [... % temperature, dry weight at end of metamorphosis (mg)
    15    30
    27    9
];
units.tempdry  = {'°C', 'mg'}; label.tempdry = {'temperature', 'dry weight'};  
bibkey.tempdry = 'Morand97';

data.tempsdur = [... % temperature, dry weight at end of metamorphosis (mg)
    15    73
    27    29.5
];
units.tempsdur  = {'°C', 'd'}; label.tempsdur = {'temperature', 'days from hatch'};  
bibkey.tempsdur = 'Morand97';

%% data from Laurila 1998
% temperature 19°C; food 1= 2% body weight; food 2 = 6% body weight
% duration of larval period counted as time from G25 to G42
% same for teh weight
data.fD = 60.8; %days from G25 to G42
units.fD = 'd'; label.fD = 'time since birth at metam';bibkey.fD = 'Laurila1998';   
temp.fD = C2K(19);  units.temp.fD = 'K'; label.temp.fD = 'temperature';
data.fW = 181.0; %days from G25 to G42
units.fW = 'd'; label.fW = 'wet weight at G42';bibkey.fW = 'Laurila1998';   
temp.fW = C2K(19);  units.temp.fW = 'K'; label.temp.fW = 'temperature';

data.fD2 = 38.2; %days from G25 to G42
units.fD2 = 'd'; label.fD2 = 'time since birth at metam';bibkey.fD2 = 'Laurila1998';   
temp.fD2 = C2K(19);  units.temp.fD2 = 'K'; label.temp.fD2 = 'temperature';
data.fW2 = 205.5; %days from G25 to G42
units.fW2 = 'd'; label.fW2 = 'wet weight at G42';bibkey.fW2 = 'Laurila1998';   
temp.fW2 = C2K(19);  units.temp.fW2 = 'K'; label.temp.fW2 = 'temperature';


%% Mikó 2017
data.tmet = 41.8; %days from G25 to G42
units.tmet = 'd'; label.tmet = 'time since birth at metam';bibkey.tmet = 'Mikó2017';   
temp.tmet = C2K(18);  units.temp.tmet = 'K'; label.temp.tmet = 'temperature';
data.Mmet = 288; %days from G25 to G42
units.Mmet = 'd'; label.Mmet = 'wet weight at G42';bibkey.Mmet = 'Mikó2017';   
temp.Mmet = C2K(18);  units.temp.Mmet = 'K'; label.temp.Mmet = 'temperature';

%% set weights for all real data
weights = setweights(data, []);

weights.tW = repelem(1,length(data.tW))';
weights.tL = repelem(1,length(data.tL))';
weights.St = repelem(6,length(data.St))';
% weights.tempsdur = repelem(1,length(data.tempsdur))';
% weights.tempdry = repelem(1,length(data.tempdry))';


%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
weights.psd.v = 3 * weights.psd.v;

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Discussion points
D1 = 'We here assume that, at metam, del_M and z change';
D2 = 'mod 2: ingestion-growth data removed, given in Jorg1988  because the time unit was not indicated and toad weight unclear';
D3 = 'mod 2: ultimate body weight is substantially increased, and sex dimorphy made explicit and Wwb reduced';
D4 = 'Males are assumed to differ from females by {p_Am} only';
metaData.discussion = struct('D1',D1, 'D2',D2, 'D3',D3, 'D4',D4);

%% Facts
F1 = ['Spends its life as a tadpole in the water and as an adult mainly on land in water rich environments. '...
      'Female toads are somewhat larger (~15 cm) than males (~13 cm). '...
      'It has a larval period of ca 58 days and it can become relatively old with a max age of ca 36 years. '... 
      'At an age of 3 - 4 years it starts reproduction. It lays a max of 8,000 eggs once a year. '...
      'During winter it hibernates and in spring massive amounts of toads can migrate to their breeding grounds. '... 
      'It has a characteristic way of laying eggs in strings.'];
metaData.bibkey.F1 = 'Wiki'; 
F2 = 'Tadpoles lose much of their weight during metamorphosis';
metaData.bibkey.F2 = 'Wiki'; 
F3 = 'length-weight males = Ww in g  = exp(-8.01 + 2.70 * 1og(SVL im mm))';
metaData.bibkey.F3 = 'ReadClar1995'; 
F4 = 'length-weight females = Ww in g  = exp(-10.1 + 3.22 * 1og(SVL im mm))';
metaData.bibkey.F4 = 'ReadClar1995'; 
metaData.facts = struct('F1',F1,'F2',F2,'F3',F3,'F4',F4);

%% Links
metaData.links.id_CoL = 'NP2D'; % Cat of Life
metaData.links.id_ITIS = '173480'; % ITIS
metaData.links.id_EoL = '333310'; % Ency of Life
metaData.links.id_Wiki = 'Bufo_bufo'; % Wikipedia
metaData.links.id_ADW = 'Bufo_bufo'; % ADW
metaData.links.id_Taxo = '47784'; % Taxonomicon
metaData.links.id_WoRMS = '1350162'; % WoRMS
metaData.links.id_amphweb = 'Bufo+bufo'; % AmphibiaWeb
metaData.links.id_AnAge = 'Bufo_bufo'; % AnAge


%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Bufo_bufo}}';
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
bibkey = 'KatzWari2003'; type = 'Article'; bib = [ ... 
'author = {Katzmann, S. and Waringer-Loschenkohl, A. and Waringer, A.}, ' ... 
'year = {2003}, ' ...
'title = {Effects of inter- and intraspecific competition on growth and development of Bufo viridis and \emph{Bufo bufo} tadpoles}, ' ...
'journal = {Limnologica}, ' ...
'volume = {33}, ' ...
'pages = {122--130}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'BrunBern2009'; type = 'Article'; bib = [ ... 
'author = {E. Brunelli and I. Bernab\`{o} and C. Bergb and K. Lundstedt-Enkel and A. Bonaccia and S. Tripepia}, ' ... 
'year = {2009}, ' ...
'title = {Environmentally relevant concentrations of endosulfan impair development, metamorphosis and behaviour in \emph{Bufo bufo} tadpoles}, ' ...
'journal = {Aquatic Toxicology}, ' ...
'volume = {91}, ' ...
'pages = {135--142}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey='Cvetkovic08'; type = 'Article'; bib = [ ... 
'author = {Cvetković, D. and Tomašević, N. and Ficetola, G. F. and Crnobrnja-Isailović, J. and Miaud, C.},' ...
"title = {Bergmann's rule in amphibians: combining demographic and ecological parameters to explain body size variation among populations in the common toad Bufo bufo},"...
"journal = {Journal of Zoological Systematics and Evolutionary Research}," ...
"volume = {47}," ...
"number = {2}," ...
"pages = {171-180}," ...
"doi = {https://doi.org/10.1111/j.1439-0469.2008.00504.x}," ...
"year = {2009}"];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey='Hemelaar88'; type = 'Article'; bib = [ ... 
'author = {Agnes Hemelaar},' ...
"title = {Age, Growth and Other Population Characteristics of Bufo bufo from Different Latitudes and Altitudes},"...
"journal = {Journal of Herpetology}," ...
"volume = {22}," ...
"number = {4}," ...
"pages = {369-388}," ...
"doi = {https://doi.org/10.2307/1564332}," ...
"year = {1988}"];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'ReadClar1995'; type = 'Article'; bib = [ ... 
'author = {C. J. Reading. R. T. Clarke}, ' ... 
'year = {1995}, ' ...
'title = {The effects of density, rainfall and environmental temperature on body condition and fecundity in the common toad, \emph{Bufo bufo}}, ' ...
'journal = {Oecologia}, ' ...
'volume = {102}, ' ...
'pages = {453-459}'];
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
%
bibkey = 'Laurila1998'; type = 'Article'; bib = [ ... 
'author = {Anssi Laurila and Jutta Kujasalo and Esa Ranta}, ' ... 
'year = {1998}, ' ...
'title = {Predator-Induced Changes in Life History in Two Anuran Tadpoles: Effects of Predator Diet}, ' ...
'journal = {Oikos}, ' ...
'volume = {83},' ...
'number = {2},' ...
'pages = {307-317}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Miko2017'; type = 'Article'; bib = [ ... 
'author = {Zsanett Mikó and János Ujszegi and Attila Hettyey}, ' ... 
'year = {2017}, ' ...
'title = {Age-dependent changes in sensitivity to a pesticide in tadpoles of the common toad (Bufo bufo)}, ' ...
'journal = {Aquatic Toxicology}, ' ...
'volume = {187}, ' ...
"doi = {https://doi.org/10.1016/j.aquatox.2017.03.016}," ...
'pages = {48-54}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'EoL'; type = 'Misc'; bib = ...
'howpublished = {\url{http://eol.org/pages/333310/overview}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'amphibiaweb'; type = 'Misc'; bib = ...
'howpublished = {\url{https://amphibiaweb.org/species/100}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];


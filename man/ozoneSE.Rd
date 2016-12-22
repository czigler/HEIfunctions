\name{ozoneSE}
\alias{ozoneSE}
\docType{data}
\title{
Ozone in the Southeast
}
\description{
Ozone in the southeast, 10/31/12
}
\usage{data(ozoneSE)}
\format{
  A data frame with 382 observations on the following 45 variables.
  \describe{
    \item{\code{Monitor_ID}}{a numeric vector}
    \item{\code{County}}{a factor with levels \code{} \code{Abbeville} \code{Ada} \code{Adair} \code{Adams} \code{Aiken} \code{Alachua} \code{Alameda} \code{Albany} \code{Albemarle} \code{Alcorn} \code{Alexander} \code{Alexandria City} \code{Allegan} \code{Allegheny} \code{Allen} \code{Amador} \code{Anchorage} \code{Anderson} \code{Andrew} \code{Androscoggin} \code{Anne Arundel} \code{Anoka} \code{Arapahoe} \code{Arlington} \code{Armstrong} \code{Aroostook} \code{Ascension} \code{Ashland} \code{Ashtabula} \code{Athens} \code{Atlantic} \code{Avery} \code{Baker} \code{Baldwin} \code{Baltimore} \code{Baltimore City} \code{Barnstable} \code{Barnwell} \code{Bay} \code{Beauregard} \code{Beaver} \code{Becker} \code{Belknap} \code{Bell} \code{Bennington} \code{Benzie} \code{Bergen} \code{Berkeley} \code{Berks} \code{Berkshire} \code{Bernalillo} \code{Berrien} \code{Bexar} \code{Bibb} \code{Billings} \code{Blair} \code{Blount} \code{Bolivar} \code{Boone} \code{Bossier} \code{Boulder} \code{Box Elder} \code{Boyd} \code{Brazoria} \code{Bremer} \code{Brevard} \code{Brewster} \code{Bristol} \code{Bronx} \code{Brookings} \code{Broward} \code{Brown} \code{Bucks} \code{Bullitt} \code{Buncombe} \code{Burke} \code{Burleigh} \code{Butler} \code{Butte} \code{Cabell} \code{Cache} \code{Caddo} \code{Calaveras} \code{Calcasieu} \code{Caldwell} \code{Callaway} \code{Calvert} \code{Cambria} \code{Camden} \code{Cameron} \code{Campbell} \code{Canadian} \code{Canyon} \code{Carbon} \code{Carlton} \code{Caroline} \code{Carroll} \code{Carson City} \code{Carter} \code{Cass} \code{Cassia} \code{Caswell} \code{Catano} \code{Cecil} \code{Cedar} \code{Centre} \code{Champaign} \code{Charles} \code{Charles City} \code{Charleston} \code{Chatham} \code{Chattooga} \code{Chautauqua} \code{Chemung} \code{Cherokee} \code{Cheshire} \code{Chester} \code{Chesterfield} \code{Chippewa} \code{Chittenden} \code{Christian} \code{Churchill} \code{Clackamas} \code{Clallam} \code{Clark} \code{Clarke} \code{Clay} \code{Clermont} \code{Cleveland} \code{Clinton} \code{Cobb} \code{Cochise} \code{Coconino} \code{Coffee} \code{Colbert} \code{Colleton} \code{Collier} \code{Collin} \code{Columbia} \code{Colusa} \code{Comanche} \code{Contra Costa} \code{Cook} \code{Coos} \code{Coweta} \code{Crittenden} \code{Crook} \code{Cumberland} \code{Custer} \code{Cuyahoga} \code{Dakota} \code{Dallas} \code{Dane} \code{Darlington} \code{Dauphin} \code{Davidson} \code{Davie} \code{Daviess} \code{Davis} \code{Dawson} \code{DeKalb} \code{Delaware} \code{Denali} \code{Denton} \code{Denver} \code{Deschutes} \code{DeSoto} \code{District of Columbia} \code{Dodge} \code{Dona Ana} \code{Door} \code{Douglas} \code{Du Page} \code{Dukes} \code{Dunn} \code{Duplin} \code{Durham} \code{Dutchess} \code{Duval} \code{Dyer} \code{East Baton Rouge} \code{Eddy} \code{Edgecombe} \code{Edgefield} \code{Edmonson} \code{Effingham} \code{El Dorado} \code{El Paso} \code{Elkhart} \code{Ellis} \code{Elmore} \code{Erie} \code{Escambia} \code{Essex} \code{Etowah} \code{Fairfax} \code{Fairfield} \code{Fauquier} \code{Fayette} \code{Flathead} \code{Florence} \code{Floyd} \code{Fond du Lac} \code{Forest} \code{Forsyth} \code{Franklin} \code{Frederick} \code{Fremont} \code{Fresno} \code{Fulton} \code{Galveston} \code{Garfield} \code{Garrett} \code{Geauga} \code{Genesee} \code{Gibson} \code{Gila} \code{Glacier} \code{Glenn} \code{Gloucester} \code{Glynn} \code{Grafton} \code{Graham} \code{Grant} \code{Granville} \code{Graves} \code{Green} \code{Greenbrier} \code{Greene} \code{Greenup} \code{Greenville} \code{Gregg} \code{Guilford} \code{Gwinnett} \code{Hamblen} \code{Hamilton} \code{Hampden} \code{Hampshire} \code{Hampton City} \code{Hancock} \code{Hanover} \code{Hardin} \code{Harford} \code{Harris} \code{Harrison} \code{Hartford} \code{Hawaii} \code{Hays} \code{Haywood} \code{Henderson} \code{Hendricks} \code{Henrico} \code{Henry} \code{Herkimer} \code{Hidalgo} \code{Highlands} \code{Hillsborough} \code{Hinds} \code{Holmes} \code{Honolulu} \code{Hood} \code{Houston} \code{Hudson} \code{Humboldt} \code{Humphreys} \code{Hunt} \code{Hunterdon} \code{Huntington} \code{Huron} \code{Iberville} \code{Imperial} \code{Indiana} \code{Ingham} \code{Inyo} \code{Jackson} \code{Jasper} \code{Jefferson} \code{Jersey} \code{Jessamine} \code{Johnson} \code{Johnston} \code{Juncos} \code{Kalamazoo} \code{Kanawha} \code{Kane} \code{Kaufman} \code{Kay} \code{Kennebec} \code{Kenosha} \code{Kent} \code{Kenton} \code{Kern} \code{Kewaunee} \code{King} \code{Kings} \code{Kleberg} \code{Klickitat} \code{Knox} \code{Kootenai} \code{Kosciusko} \code{La Crosse} \code{La Paz} \code{La Plata} \code{La Porte} \code{Lackawanna} \code{Lafayette} \code{Lafourche} \code{Lake} \code{Lancaster} \code{Lane} \code{Larimer} \code{Larue} \code{Lauderdale} \code{Lawrence} \code{Lea} \code{Leavenworth} \code{Lee} \code{Leelanau} \code{Lehigh} \code{Lenawee} \code{Lenoir} \code{Leon} \code{Lewis} \code{Licking} \code{Limestone} \code{Lincoln} \code{Linn} \code{Litchfield} \code{Livingston} \code{Lorain} \code{Los Angeles} \code{Loudon} \code{Loudoun} \code{Love} \code{Lucas} \code{Luna} \code{Luzerne} \code{Lycoming} \code{Lyon} \code{Macomb} \code{Macon} \code{Macoupin} \code{Madera} \code{Madison} \code{Mahoning} \code{Manatee} \code{Manitowoc} \code{Marathon} \code{Maricopa} \code{Marin} \code{Marion} \code{Mariposa} \code{Marshall} \code{Martin} \code{Mason} \code{McCracken} \code{McCurtain} \code{McHenry} \code{McKenzie} \code{McLean} \code{McLennan} \code{Meade} \code{Mecklenburg} \code{Medina} \code{Meigs} \code{Mendocino} \code{Merced} \code{Mercer} \code{Merrimack} \code{Mesa} \code{Miami} \code{Miami-Dade} \code{Middlesex} \code{Mille Lacs} \code{Milwaukee} \code{Minnehaha} \code{Missaukee} \code{Missoula} \code{Mobile} \code{Monmouth} \code{Mono} \code{Monongalia} \code{Monroe} \code{Monterey} \code{Montezuma} \code{Montgomery} \code{Morgan} \code{Morris} \code{Muhlenberg} \code{Multnomah} \code{Murray} \code{Muscogee} \code{Muskegon} \code{Napa} \code{Navajo} \code{Nevada} \code{New Castle} \code{New Hanover} \code{New Haven} \code{New London} \code{New York} \code{Newport News City} \code{Newton} \code{Niagara} \code{Norfolk} \code{Northampton} \code{Nueces} \code{Oakland} \code{Obion} \code{Ocean} \code{Oconee} \code{Ohio} \code{Okaloosa} \code{Oklahoma} \code{Oldham} \code{Oliver} \code{Olmsted} \code{Oneida} \code{Onondaga} \code{Orange} \code{Orangeburg} \code{Orleans} \code{Osceola} \code{Oswego} \code{Ottawa} \code{Ouachita} \code{Outagamie} \code{Oxford} \code{Ozaukee} \code{Page} \code{Palm Beach} \code{Palo Alto} \code{Parker} \code{Pasco} \code{Passaic} \code{Paulding} \code{Pennington} \code{Penobscot} \code{Peoria} \code{Perry} \code{Person} \code{Philadelphia} \code{Pickens} \code{Pierce} \code{Pike} \code{Pima} \code{Pinal} \code{Pinellas} \code{Piscataquis} \code{Pitkin} \code{Pitt} \code{Pittsburg} \code{Placer} \code{Platte} \code{Plumas} \code{Pointe Coupee} \code{Polk} \code{Portage} \code{Porter} \code{Posey} \code{Powder River} \code{Preble} \code{Prince Georges} \code{Prince William} \code{Providence} \code{Pulaski} \code{Putnam} \code{Queens} \code{Racine} \code{Randolph} \code{Rensselaer} \code{Richland} \code{Richmond} \code{Rio Blanco} \code{Riverside} \code{Roanoke} \code{Rock} \code{Rock Island} \code{Rockbridge} \code{Rockdale} \code{Rockingham} \code{Rockland} \code{Rockwall} \code{Rosebud} \code{Rowan} \code{Russell} \code{Rutherford} \code{Sacramento} \code{Sagadahoc} \code{Saint Bernard} \code{Saint Charles} \code{Saint Clair} \code{Saint Croix} \code{Saint James} \code{Saint John} \code{Saint Louis} \code{Saint Louis City} \code{Saint Lucie} \code{Saint Mary} \code{Saint Tammany} \code{Sainte Genevieve} \code{Salt Lake} \code{San Benito} \code{San Bernardino} \code{San Diego} \code{San Francisco} \code{San Joaquin} \code{San Juan} \code{San Luis Obispo} \code{San Mateo} \code{Sandoval} \code{Sangamon} \code{Santa Barbara} \code{Santa Clara} \code{Santa Cruz} \code{Santa Fe} \code{Santa Rosa} \code{Sarasota} \code{Saratoga} \code{Sauk} \code{Schenectady} \code{Schoolcraft} \code{Scott} \code{Scotts Bluff} \code{Sedgwick} \code{Seminole} \code{Sevier} \code{Shasta} \code{Shawnee} \code{Sheboygan} \code{Shelby} \code{Simpson} \code{Sioux} \code{Siskiyou} \code{Skagit} \code{Smith} \code{Solano} \code{Sonoma} \code{Spartanburg} \code{Spokane} \code{St John the Baptist} \code{St Joseph} \code{Stafford} \code{Stanislaus} \code{Stark} \code{Stearns} \code{Steuben} \code{Story} \code{Strafford} \code{Sublette} \code{Suffolk} \code{Suffolk City} \code{Sullivan} \code{Summit} \code{Sumner} \code{Sumter} \code{Sussex} \code{Sutter} \code{Swain} \code{Sweetwater} \code{Talladega} \code{Taney} \code{Tarrant} \code{Tehama} \code{Teton} \code{Thurston} \code{Tolland} \code{Tooele} \code{Travis} \code{Trego} \code{Trigg} \code{Trumbull} \code{Tulare} \code{Tulsa} \code{Tuolumne} \code{Tuscaloosa} \code{Uinta} \code{Uintah} \code{Ulster} \code{Umatilla} \code{Union} \code{Utah} \code{Valencia} \code{Van Buren} \code{Vanderburgh} \code{Ventura} \code{Vernon} \code{Victoria} \code{Vigo} \code{Vilas} \code{Volusia} \code{Wake} \code{Wakulla} \code{Walker} \code{Walworth} \code{Warren} \code{Warrick} \code{Washington} \code{Washoe} \code{Washtenaw} \code{Waukesha} \code{Wayne} \code{Webb} \code{Weber} \code{Webster} \code{Weld} \code{West Baton Rouge} \code{Westchester} \code{Westmoreland} \code{Whatcom} \code{White Pine} \code{Will} \code{Williamsburg} \code{Williamson} \code{Wilson} \code{Winnebago} \code{Wood} \code{Worcester} \code{Wright} \code{Wyandotte} \code{Wythe} \code{Yancey} \code{Yavapai} \code{Yellowstone} \code{Yolo} \code{York} \code{Yuma}}
    \item{\code{State}}{a factor with levels \code{} \code{AK} \code{AL} \code{AR} \code{AZ} \code{CA} \code{CO} \code{CT} \code{DC} \code{DE} \code{FL} \code{GA} \code{HI} \code{IA} \code{ID} \code{IL} \code{IN} \code{KS} \code{KY} \code{LA} \code{MA} \code{MD} \code{ME} \code{MI} \code{MN} \code{MO} \code{MS} \code{MT} \code{NC} \code{ND} \code{NE} \code{NH} \code{NJ} \code{NM} \code{NV} \code{NY} \code{OH} \code{OK} \code{OR} \code{PA} \code{PR} \code{RI} \code{SC} \code{SD} \code{TN} \code{TX} \code{UT} \code{VA} \code{VI} \code{VT} \code{WA} \code{WI} \code{WV} \code{WY}}
    \item{\code{Longitude}}{a numeric vector}
    \item{\code{Latitude}}{a numeric vector}
    \item{\code{location}}{a factor with levels \code{RURAL} \code{SUBURBAN} \code{UNKNOWN} \code{URBAN AND CENTER CITY}}
    \item{\code{a}}{a numeric vector}
    \item{\code{Median_income}}{a numeric vector}
    \item{\code{HS_rate}}{a numeric vector}
    \item{\code{Urban_rate}}{a numeric vector}
    \item{\code{Migration_5_year_rate}}{a numeric vector}
    \item{\code{Current_Smoking_rate_weighted}}{a numeric vector}
    \item{\code{HOUSEUNIT}}{a numeric vector}
    \item{\code{POPULATION}}{a numeric vector}
    \item{\code{Age_med2004}}{a numeric vector}
    \item{\code{Female_med2004}}{a numeric vector}
    \item{\code{White_med2004}}{a numeric vector}
    \item{\code{Black_med2004}}{a numeric vector}
    \item{\code{allcause_death2004}}{a numeric vector}
    \item{\code{Person_years2004}}{a numeric vector}
    \item{\code{Total_den2004}}{a numeric vector}
    \item{\code{Resp_adm2004}}{a numeric vector}
    \item{\code{CVD_adm2004}}{a numeric vector}
    \item{\code{allcause_death2010}}{a numeric vector}
    \item{\code{Person_years2010}}{a numeric vector}
    \item{\code{Total_den2010}}{a numeric vector}
    \item{\code{Resp_adm2010}}{a numeric vector}
    \item{\code{CVD_adm2010}}{a numeric vector}
    \item{\code{o38hrmaxbase}}{a numeric vector}
    \item{\code{o38hrmaxfu}}{a numeric vector}
    \item{\code{o38hr4thbase}}{a numeric vector}
    \item{\code{o38hr4thfu}}{a numeric vector}
    \item{\code{o324hrbase}}{a numeric vector}
    \item{\code{o324hrfu}}{a numeric vector}
    \item{\code{weight_cdc}}{a numeric vector}
    \item{\code{smoking_cdc}}{a numeric vector}
    \item{\code{adulthouse_cdc}}{a numeric vector}
    \item{\code{exercise_cdc}}{a numeric vector}
    \item{\code{asthma_cdc}}{a numeric vector}
    \item{\code{HTN_cdc}}{a numeric vector}
    \item{\code{diabetes_cdc}}{a numeric vector}
    \item{\code{BMIover_cdc}}{a numeric vector}
    \item{\code{BMIob_cdc}}{a numeric vector}
    \item{\code{genhealth_cdc}}{a numeric vector}
    \item{\code{nodoc_cdc}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ozoneSE)
## maybe str(ozoneSE) ; plot(ozoneSE) ...
}
\keyword{datasets}

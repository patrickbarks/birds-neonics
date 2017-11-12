birds-neonics
==============

Data and code for analyses described in the manuscript "Declines in insectivorous birds in the United States not associated with the use of neonicotinoid insecticides", by Patrick Barks and Paul Galpern.

Some of the files used in our analysis are not included in this repository due to size limitations. These include:
1. `data/bbs-allstates-2015.txt`
2. `data/epest-merge.csv`
3. `data/nlcd_2011_landcover_2011_edition_2014_10_10/`
4. `stanfit/`

File #1 (`bbs-allstates-2015.txt`) is a concatenation of Breeding Bird Survey data available [here](ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/Archivefiles/Version2015v1/States/), in separate tables for each state/province. It can be downloaded via Google Drive link, or created using the following shell commmands on Mac OS X:

```
cd data/bbs_raw_2015/States/
ftp -i ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/Archivefiles/Version2015v1/States/
mget *.zip
exit
unzip "*.zip"
rm *.zip
head -1 Alabama.csv > ../../bbs-allstates-2015.txt
tail -n +2 -q *.csv >> ../../bbs-allstates-2015.txt
```

File #2 can be produced using `r/build-epest.R`. Folder #3 can be downloaded [here](https://www.mrlc.gov/nlcd06_data.php). The files within folder #4 are model objects totalling 39 GB, that can be downloaded via Google Drive link. (Note our main analyses (`r/model-aggregate.R` and `r/plot.R`) can still be run without the model objects).


### Sources for raw data files (folder `data/`)

##### North American Breeding Bird Survey (2015 version)
[ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/Archivefiles/Version2015v1/](ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/Archivefiles/Version2015v1/)
- `bbs_raw_2015/States/`
- `bbs_raw_2015/BCR.csv`
- `bbs_raw_2015/RegionCodes.csv`
- `bbs_raw_2015/routes.csv` (within Routes.zip)
- `bbs_raw_2015/weather.csv` (within Weather.zip)

[https://www.pwrc.usgs.gov/bbs/stratanames/index.html](https://www.pwrc.usgs.gov/bbs/stratanames/index.html)
- `bbs_raw_2015/bbs_stata_info.csv` (csv file derived manually from info on webpage)

[https://www.mbr-pwrc.usgs.gov/bbs/geographic_information/geographic_information_products_.htm](https://www.mbr-pwrc.usgs.gov/bbs/geographic_information/geographic_information_products_.htm)
- `bbs_raw_2015/bbsrte_2012_alb/` (note: apparently no longer available)

[https://www.pwrc.usgs.gov/bba/index.cfm?fa=bba.getData](https://www.pwrc.usgs.gov/bba/index.cfm?fa=bba.getData)
- `bbs_raw_2015/bcr_shp/` (within bcr_shp.zip)

##### EPEST High
[https://pubs.usgs.gov/ds/752/](https://pubs.usgs.gov/ds/752/)
- `epest_raw/1992_2009/` (see 14 text files on right side of page)

[https://pubs.usgs.gov/ds/0907/](https://pubs.usgs.gov/ds/0907/)
- `epest_raw/2008_2012/` (see 14 text files on right side of page)

##### National Land Cover Database (2011)
[https://www.mrlc.gov/nlcd06_data.php](https://www.mrlc.gov/nlcd06_data.php)
- `nlcd_2011_landcover_2011_edition_2014_10_10/`

##### American counties shapefile
[https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html](https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html)
- `gz_2010_us_050_00_500k/` (path `2010 Census` > gz_2010_us_050_00_500k.zip)

##### Neonicotinoid use estimates for the Netherlands
[http://statline.cbs.nl/statweb](http://statline.cbs.nl/statweb)
- `neonic-use-netherlands.csv` (csv derived manually via path: `Theme` > `Agriculture` > `Pesticides` > `Crop protection; active ingredient` > `Amount of use` > `Use per year` > [select compounds and periods of interest])


### Summary of R scripts (folder `r/`)

##### `build-bbs.R`
Arrange the various Breeding Bird Survey data tables into a single data frame (saved as `data/bbs-use.RData`), which is subset to our 28 species of interest, and which includes counts of zero (i.e. instances where a given species was not observed during a particular survey).

##### `cartography.R`
Geographic analyses to assign BBS routes to counties (`data/route-county-start.csv`), and estimate the proportion of land area devoted to cultivated crops around each BBS route (`data/landcover-route.csv`). Also creates a grid of quadrats laid overtop the contiguous USA, and assigns BBS routes to quadrats (`data/route-quadrat.csv`). Also produces various shapefiles for plotting.

##### `build-epest.R`
Arrange the various EPEST data tables into a single data frame (saved as `data/epest-merge.csv`). Also produces data frames for total neonicotinoid use in each county and year scaled by county area (`data/neonic-county-year.csv`), average annual neonicotinoid use in each county over the period 2005-2012 scaled by county area (`data/neonic-county-mean-2005-2012.csv`), and average annual imidacloprid use in each county over the period 2005-2012 scaled by county area (`data/imidacloprid-county-mean-2005-2012.csv`).

##### `model-spp-[TYPE].R`
Fits species-level models using RStan library, and saves resulting stanfit objects to the folder `stanfit/`. These scripts were run on Amazon Web Service EC2 instances (type c4.xlarge). For a given model type and time period, the run-time for the full list of 28 species on a single instance was on the order of 72 hours.

##### `get-post-summary-[TYPE].R`
Extracts model diagnostics and posterior samples for parameters of interest, from each of the species-level models (i.e. from the stanfit objects in the folder `stanfit/`). Resulting summaries are save as `analysis/post-summary-spp-[TYPE].csv`.

##### `model-aggregate.R`
Fits aggregation models (i.e. models to aggregate parameters of interest across species and regions) using RStan, and saves resulting summaries as `analysis/post-summary-agg-[TYPE].RData`.

##### `plot.R`
Generates all figures within manuscript.

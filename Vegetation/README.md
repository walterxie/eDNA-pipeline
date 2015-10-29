# Little Barrier Island Project

In development

## Files

- Raw data: [LBI\_Full\_plant\_data\_v1.7.xls](https://raw.githubusercontent.com/walterxie/eDNA-pipeline/hauturu/Vegetation/LBI_Full_plant_data_v1.7.xls)
  
- Vegetation plot metadata: [LBI\_veg\_metadata.txt](https://raw.githubusercontent.com/walterxie/eDNA-pipeline/hauturu/Vegetation/LBI_veg_metadata.txt)
  
- Script to create SBA matrix from DBH: [DBHtoSBA.r](https://raw.githubusercontent.com/walterxie/eDNA-pipeline/hauturu/Vegetation/DBHtoSBA.r)
  
- DBH: LBI\_Tree\_DBH\_Values\_All\_Trees_Saplings.csv

Column names should not have any white space. Example file format:
```
PlotNo	PlotName	Species	StemCount	DBH
1	Plot01	BRAREP	1	1.3
1	Plot01	BRAREP	1	4
```

- SBA (cm<sup>2</sup>/m<sup>2</sup>) matrix: 
[LBI\_Trees\_Saplings\_SBA.csv](https://raw.githubusercontent.com/walterxie/eDNA-pipeline/hauturu/Vegetation/LBI_Trees_Saplings_SBA.csv)


## Vegetation Plot Metadata 

* Plot No.: Designated number of the vegetation plot.				
* Plot Name: Designated name of the vegetation plot.				
* Date Inventoried: Date the plot was inventoried, or start date for those plots that required two consecutive days.				
* Elevation: Elevation of the vegetation plots 'P' corner, in metres.				
* Slope: Slope of the plot, in degrees.				
* Aspect: Aspect of the plot, in degrees.				
* Physiography: Plot physiography, based on broad categories.				
* Drainage: Drainage of the plots soil - ranked poor, moderate, or good.				
* Latitude: Decimal degree based latitude of the plots 'P' corner.				
* Longitude: Decimal degree based longitude of the plots 'P' corner.				
* Lat-Long Datum: Geodetic datum associated with the latitudes and longitudes.				
* Easting: Easting of the plots 'P' corner, using the New Zealand Map Grid (NZMG) projection.				
* Northing: Northing of the plots 'P' corner, using the New Zealand Map Grid (NZMG) projection.				
* East-North Datum: Geodetic datum associated with the eastings and northings.				


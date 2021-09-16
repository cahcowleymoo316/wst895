drop table #tmp
SELECT a.[Suburb_Id]
, b.suburb
, b.municipality
, cast(a.[Crime_Res_1000hh] as numeric) as crime_res_1000hh_use
, cast(a.[Crime_Vehicle_1000hh] as numeric) as crime_vehicle_1000hh_use
, c.Adults
, b.x
, b.y
, b.SP_GEOMETRY.STAsBinary() as Sub_WKB
into #tmp
FROM [BAS_Insurance].[dbo].[WS_Suburb_Metrics] a
inner join spatial_base..ls_suburbs b
on a.suburb_id = b.suburb_id
inner join BAS_Insurance..suburb_demographics c
on a.suburb_id = c.suburb_id
where b.PROVINCE= 'GAUTENG' 
and a.[Crime_Res_1000hh] <> '<4'
and a.[Crime_Res_1000hh] is not null
and a.crime_vehicle_1000hh <> '<4'
and a.Crime_Vehicle_1000hh <>'>500'
and a.crime_vehicle_1000hh is not null
order by b.suburb


select * into usr_charl..gauteng_suburbs from #tmp where crime_vehicle_1000hh_use < 60 and crime_res_1000hh_use < 60


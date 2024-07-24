
-- 26 - TOTAL EGGS LAID 
SELECT SUM(Total_Eggs_Laid) AS Total_Eggs_Laid FROM worm_summary

-- 26m TOTAL EGGS LAID Mass
SELECT SUM(Total_Eggs_Laid) * 65 AS Total_Eggs_Laid_MASS FROM worm_summary

-- 28 - Count number of eggs removed by culling
SELECT COUNT(*) from worms WHERE stage = "dead" AND Notes LIKE "%cull%" AND Mass <= 65


-- 28m Mass of Eggs removed by culling
SELECT COUNT(*) * 65 from worms WHERE stage = "dead" AND Notes LIKE "%cull%" AND Mass <= 65

-- 30 Count number of eggs hatched
SELECT COUNT (*) as Hatched FROM worm_summary WHERE Larva_span_days is NOT NULL


-- 30 Mass of eggs hatched
SELECT COUNT (*) * 65 as Hatched FROM worm_summary WHERE Larva_span_days is NOT NULL

-- 31 - Count eggs into L1 Arrest
-- TODO Not yet implemented --


-- V Worms all

-- 14 - Count unique number of worms
SELECT COUNT(DISTINCT Worm_Name) AS unique_worm_count
FROM worm_summary;


-- 14a - Count total body mass of worms in ng
SELECT SUM(Total_Body_Mass) AS total_body_mass_sum
FROM worm_summary


--17 - Count born as dauer by counting the number of worms that have a Dauer_span_days but not a Larva_span_days
SELECT COUNT(*) AS BORN_DAUER
FROM worm_summary
WHERE Dauer_span_days IS NOT NULL and Larva_span_days IS NULL

-- 17a Count by mass of born as dauer
SELECT SUM(Total_Body_Mass) AS BORN_DAUER_MASS_NG
FROM worm_summary
WHERE Dauer_span_days IS NOT NULL and Larva_span_days IS NULL


-- 18 - Count born as egg by counting the number of worms that have a Larva_span_days
SELECT SUM(Total_Body_Mass) AS BORN_EGG
FROM worm_summary
WHERE  Larva_span_days IS NOT NULL

-- 18a - Count mass of worms born as egg by counting the number of worms that have a Larva_span_days
SELECT COUNT(*) AS BORN_EGG_MASS
FROM worm_summary
WHERE  Larva_span_days IS NOT NULL


-- 19  Count all dead worms:
SELECT COUNT(*) AS DEAD_WORM_COUNT
FROM worms
WHERE Notes LIKE "%Cause of Death%"

-- 19a Count all worms died by culled
SELECT COUNT(*) AS DEAD_WORM_CULLED_COUNT
FROM worms
WHERE Notes LIKE "%culled%"

-- 19a_m Count all worms died by culled mass
SELECT SUM(Mass) AS DEAD_WORM_CULLED_MASS
FROM worms
WHERE Notes LIKE "%culled%"

-- 19b Count all worms died by starvation
SELECT COUNT(*) AS DEAD_WORM_STARVE_COUNT
FROM worms
WHERE Notes LIKE "%starve%"

-- 19b_m Count all worms died by starvation mass
SELECT SUM(Mass) AS DEAD_WORM_STARVE_MASS
FROM worms
WHERE Notes LIKE "%starve%"


-- 19c Count all worms died by old age
SELECT COUNT(*) AS DEAD_WORM_COUNT
FROM worms
WHERE Notes LIKE "%old age%"

-- 19c_m Count all worms died by old age
SELECT SUM(Mass) AS DEAD_WORM_MASS
FROM worms
WHERE Notes LIKE "%old age%"


-- 19 Count all dead mass
SELECT SUM(Mass) AS DEAD_MASS
FROM worms WHERE Notes LIKE "%Cause of Death%"

-- 20 - Count all worms that are alive
------- NOT IMPLEMENTED IN SQL ------


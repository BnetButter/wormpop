SELECT (SELECT SUM(larva_to_adult) from stage_transition) as 'Adult  In',
	(SELECT SUM(adult_culled_ind) from dynamic_stage_transition) as 'Adult Culled',
	(SELECT SUM(adult_to_bag) from stage_transition) as 'Adult to Parlad',
	(SELECT COUNT(*) FROM worms WHERE Notes like '%old%') as 'Adults Reached Old Age',
	(SELECT Adult_end_of_simulation_ind FROM dynamic_stage_transition ORDER BY timestep DESC LIMIT 1) as 'Adults Remaining',
	(SELECT COUNT(*) FROM worm_summary WHERE Total_Eggs_Laid > 0) as 'Reproductive Active Adults',
	(SELECT COUNT(*) FROM worm_summary WHERE Total_Eggs_Laid > 0 AND Cause_of_Death = 'culled') as 'Reproductively Active Adults Died of Culled',
	(SELECT COUNT(*) FROM worm_summary WHERE Total_Eggs_Laid > 0  AND Cause_of_Death = 'old age') as 'Reproductively Active Adults Died of Old Age'
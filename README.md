# Project Lumen — The Great Backyard Firefly Count

## The story

Fireflies are vanishing from the suburbs, and the Ridgeway Nature Trust wants to know whether
anything can be done about it. In June 2025 the Trust recruited **25 volunteer households** and
asked each of them to do the same thing for **20 consecutive nights**: walk into the backyard
after dusk, start a timer, and count every firefly flash they could see.

Volunteers were given a weather station and a soil probe, and logged four readings per night
before they started counting: air temperature, relative humidity, a soil-moisture index from the
probe, and how many minutes they actually stayed outside. (Some nights people got eaten alive by
mosquitoes and quit after ten minutes. Some nights they stayed out for an hour.)

Halfway through the season the Trust ran its flagship intervention: a **free native-plant kit** —
switchgrass, goldenrod, leaf-litter mulch — delivered to every household, with a volunteer crew to
plant it. The crews worked through the neighbourhood over several weeks, so different yards were
converted on different nights. `native_planting_pct` records the percentage of each yard planted
with natives on each survey night, so it steps up the evening a yard was converted.

The Trust's board meets in the autumn and wants an answer to one question:

> **Do native plantings bring the fireflies back? And what else drives how many fireflies a yard sees?**

There is a sub-plot. The Trust's grant officer has already seen a summary table someone made in a
spreadsheet, and is telling everyone the native-plant programme was a disaster. The board would
like a second opinion.

## The data

`project_lumen_firefly_counts.csv` — 500 rows, one row per yard-night (25 yards × 20 nights).

| Column | Type | Description |
|---|---|---|
| `yard_id` | string | Household identifier, `Y01`–`Y25`. Each yard appears on all 20 nights. |
| `survey_date` | date | Calendar date of the survey. |
| `survey_night` | integer | Night of the study, 1–20. |
| `temp_c` | numeric | Air temperature at the start of the count, °C. |
| `humidity_pct` | numeric | Relative humidity at the start of the count, %. |
| `soil_moisture_index` | numeric | Soil-probe reading, 0 (bone dry) to 1 (saturated). |
| `native_planting_pct` | numeric | Percentage of the yard planted with native species on that night, 0–100. |
| `survey_minutes` | integer | Minutes the volunteer spent counting. |
| `fireflies_counted` | integer | **Outcome.** Number of firefly flashes counted. |

Five predictors, one outcome, plus the yard and date labels from the study log.

## Notes from the field coordinator

- Nobody missed a night; the panel is complete and balanced.
- The weather station and soil probe are the same model in every yard, calibrated before the season.
- Volunteers were not told what the Trust expected to find.
- Counting effort was never standardised — that was a known compromise of using volunteers.

Good luck. The board is expecting a figure.

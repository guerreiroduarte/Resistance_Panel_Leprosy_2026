Read depth and integrity rate
================
2026-07-16

### Read depth

Comparing the total read depth and the full-length read depth will show
us the sequencing integrity.

``` r
depth <- read.csv(here("data/processed_data/amplicon_stats.csv"),
                  check.names = F) %>%
    select(SAMPLE, GENE, FVDEPTH, FDEPTH)

depth %>%
    mutate(integrity_rate = FVDEPTH / FDEPTH) %>%
    group_by(GENE) %>%
    summarise(
        mean_depth = mean(FDEPTH),
        min_depth = min(FDEPTH),
        max_depth = max(FDEPTH),
        mean_flength = mean(FVDEPTH),
        min_flength = min(FVDEPTH),
        max_flength = max(FVDEPTH),
        mean_integrity = mean(integrity_rate) * 100,
        sd_integrity = sd(integrity_rate) * 100
    ) %>%
    ungroup() %>%
    mutate(
        depth_label = sprintf("%.1f (%.1f - %.1f)", mean_depth, min_depth, max_depth),
        full_label = sprintf("%.1f (%.1f - %.1f)", mean_flength, min_flength, max_flength),
        integrity_label = sprintf("%.1f%% ± %.1f%%", mean_integrity, sd_integrity)
    ) %>%
    select(GENE, depth_label, full_label, integrity_label) %>% 
    rename(
        `Amplicon target` = GENE,
        `Mean read depth` = depth_label,
        `Mean read depth (full-length)` = full_label,
        `Mean integrity rate` = integrity_label
    ) %>% 
    kable() %>% 
    kable_classic() %>% 
    column_spec(1, italic = T)
```

<table class=" lightable-classic" style="color: black; font-family: &quot;Arial Narrow&quot;, &quot;Source Sans Pro&quot;, sans-serif; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Amplicon target
</th>

<th style="text-align:left;">

Mean read depth
</th>

<th style="text-align:left;">

Mean read depth (full-length)
</th>

<th style="text-align:left;">

Mean integrity rate
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;font-style: italic;">

23S rRNA I
</td>

<td style="text-align:left;">

5364.0 (826.9 - 26209.0)
</td>

<td style="text-align:left;">

5225.1 (785.0 - 25974.0)
</td>

<td style="text-align:left;">

95.6% ± 2.5%
</td>

</tr>

<tr>

<td style="text-align:left;font-style: italic;">

23S rRNA II
</td>

<td style="text-align:left;">

2223.7 (452.2 - 7961.6)
</td>

<td style="text-align:left;">

2076.6 (414.0 - 7702.0)
</td>

<td style="text-align:left;">

90.2% ± 4.9%
</td>

</tr>

<tr>

<td style="text-align:left;font-style: italic;">

folP1
</td>

<td style="text-align:left;">

2667.0 (20.4 - 12910.6)
</td>

<td style="text-align:left;">

2476.4 (18.0 - 12512.0)
</td>

<td style="text-align:left;">

87.2% ± 10.2%
</td>

</tr>

<tr>

<td style="text-align:left;font-style: italic;">

folP2
</td>

<td style="text-align:left;">

1332.0 (7.8 - 7918.9)
</td>

<td style="text-align:left;">

1196.4 (6.0 - 7443.0)
</td>

<td style="text-align:left;">

80.5% ± 7.8%
</td>

</tr>

<tr>

<td style="text-align:left;font-style: italic;">

gyrA
</td>

<td style="text-align:left;">

6951.6 (109.7 - 24007.5)
</td>

<td style="text-align:left;">

6495.7 (105.0 - 23215.0)
</td>

<td style="text-align:left;">

92.0% ± 3.5%
</td>

</tr>

<tr>

<td style="text-align:left;font-style: italic;">

gyrB
</td>

<td style="text-align:left;">

5058.2 (215.8 - 23330.6)
</td>

<td style="text-align:left;">

4835.8 (186.0 - 22894.0)
</td>

<td style="text-align:left;">

93.1% ± 3.8%
</td>

</tr>

<tr>

<td style="text-align:left;font-style: italic;">

rpoB
</td>

<td style="text-align:left;">

1365.2 (298.0 - 4937.3)
</td>

<td style="text-align:left;">

1284.2 (273.0 - 4679.0)
</td>

<td style="text-align:left;">

91.7% ± 4.1%
</td>

</tr>

</tbody>

</table>

Overall, all amplicons generated enough reads for variant analysis.
However, the folate pathway-related genes (*folP1* and *folP2*) could
not be efficiently sequenced in all samples – which is visible by their
minimum depth across our coort.

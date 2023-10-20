# Seurat Object Maker

## Screenshots of current state

Selecting directory, with data validity checks:

![](screenshots/page1.png)

Filtering page, setting thresholds with sliders. Violin plots and nbr cells left after filtering are updated in realtime.

![](screenshots/page2_a.png)

The thresholds limits are updated in real-time, press "Confirm filtering thresholds" to continue:

![](screenshots/page2_b.png)

Showing a loader while processing the data..

![](screenshots/page2_processing.png)

Results page:

![](screenshots/page3_a.png)

More of results, differentially expressed genes per cluster (an accordion, can open for a cluster of interest at a time (to not swamp page if many clusters)):
![](screenshots/page3_b.png)

Having pressed the "Download Seurat Object" button, the user is prompted to save the object as a `.rds` file (a filename including the date is suggested by default):

![](screenshots/page3_c.png)

---

Note that "debugging" values are still printed etc, the webapp is not mature yet. Code structure is also "prototype-y". But the app is functional.
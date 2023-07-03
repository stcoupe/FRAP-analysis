
<h1 align="center">FRAP-analysis</h1>

  <p align="center">
    This is a simple set of functions for analyzing FRAP time-course experiments performed on biomolecular condensates/complex coacervates
    <br />
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">View Demo</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>



<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

Requires Python (>= 3.7.0)

Python package dependencies
* scikit-image (>= 0.18)
* matplotlib
* numpy 
* scipy (Only required for curve fitting within the demo)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Usage

At the moment, this is a set a functions which can be imported and run on a set of FRAP timecourse images.

[FUTURE] : Make this easily run from the command line.

Three basic functions are present in the generate_frap_curve_functions:
* ```generate_frap_curves``` - This is the only function users need to directly interact with unless troubleshooting burn spot identification
* ```get_burn_center```      - Called within generate_frap_curves, and uses the images before and after photobleaching to identify the burn spot
* ```mean_circle```          - Called within get_burn_center. Gets average of all pixels within some radius of a specified point

To run the analysis on your own data:
1) Crop around each droplet which has been FRAPped. These scripts were built around having multiple droplets from a single field of view which were FRAPped as part of the same timecourse. As such, the time-keeping (time interval and number of frames) are assumed the same across all droplets being analyzed simultaneously.
2) Save individual droplets to a folder within the intended output directory called "signal_rois" such that the file structure is '/save_dir/signal_rois/xx.tif', where xx is the filename for an individual FRAP timecourse tiffstack.
3) Import 'generate_frap_curves':

   ```python
   from generate_frap_curve_fxns import generate_frap_curves
   ```
4) Specify ```read_directory```, ```save_directory```, ```time_array```, ```frap_frame```, ```frap_radius```, and ```save_output```, where
   * ```read_directory``` : directory the individual droplet FRAP data is stored. See Step 2. Nothing but the droplet files should be in this folder
   * ```save_directory``` : directory things will be saved in
   * ```time_array``` : time keeping for each frame of the timecourse. A time of 0 should correspond to the time the first FRAP event occurs (for fitting purposes)
   * ```frap_frame``` : frame the FRAP event first appears
   * ```frap_radius``` : radius of the FRAP spot in pixels
   * ```save_output``` : option to save images of the identified FRAP spot (1 is yes)
5) Run analysis. The output ```normalized_intensities``` will be an n-by-t array, with n-droplets rows and t-timepoints columns
   ```python
   normalized_intensities = generate_frap_curves(read_directory, save_directory, time_array, frap_frame, frap_radius, save_output)
   ```

7) Plot traces and fit curves to model of choice. See demo.

[FUTURE] : Include library of analysis scripts for curve fitting which incorporate different diffusion models.


<p align="right">(<a href="#readme-top">back to top</a>)</p>





<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

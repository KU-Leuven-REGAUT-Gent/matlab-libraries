# MATLAB Libraries

Version 2.0

[![Build Status](https://travis-ci.org/fredericdepuydt/matlab-libraries.svg?branch=master)](https://travis-ci.org/fredericdepuydt/matlab-libraries) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
***
***
***
*** To avoid retyping too much info. Do a search and replace for the following:
*** github_username, repo_name, twitter_handle, email, project_title, project_description
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

<br />
<p align="center">
  <h3 align="center">MATLAB libraries</h3>


  <p align="center">
    KU Leuven Gent Automation MATLAB libraries for scope files and PROFINET decoding, etc.
    <br />
    <a href="https://github.com/DimitriDeSchuyter/matlab-libraries"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/DimitriDeSchuyter/matlab-libraries/issues">Report Bug</a>
    ·
    <a href="https://github.com/DimitriDeSchuyter/matlab-libraries/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-repository">About The Repository</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
     <li>
      <a href="#created-classes-and-functions">Created Classes and functions</a>
      <ul>
        <li><a href="#scope-class">Scope class</a></li>
        <li><a href="#eth-class">Eth class</a></li>
      </ul>
    </li>
    <li><a href="#license">License</a></li>
  </ol>
</details>



<!-- ABOUT THE REPOSITORY -->
## About The repository
This repository contains matlab script to read Tektronix DPO4054B,MSO5054B and DPO2024 scopes, PROFINET decoding and a own plot function.


### Built With

* MATLAB R2020a 



<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/DimitriDeSchuyter/matlab-libraries.git
   ```
2. Add the folder and underlaying folders to the MATLAB path (Home--> Environment--Set Path)


<!-- CREATED CLASSES AND FUNCTIONS -->
## Created classes and functions
### Scope class
All the needed files are inside the scope folder.

Commands to create a new scope object for the following scope types:
* DPO4054B and DPO2024: scope.isfread(_path to file_,_verbose mode_)
* MSO5054B: scope.wfmread(_path to file_,_verbose mode_)

#### Properties of the class
This will create a scope object with the following properties inside:

<table>
<tr>
<th> Property name </th>
<th> Implemented </th>
</tr>
  <tr><td>model</td><td>Not yet</td></tr> 
  <tr><td>fileName</td><td>Yes</td></tr> 
  <tr><td>firmware_version</td><td>Not yet</td></tr> 
  <tr><td>waveform_type</td><td>Yes</td></tr> 
  <tr><td>point_format</td><td>Not yet</td></tr> 
  <tr><td>horizontal_units</td><td>Yes</td></tr> 
  <tr><td>horizontal_scale</td><td>Not yet</td></tr> 
  <tr><td>horizontal_delay</td><td>Not yet</td></tr> 
  <tr><td>sample_interval</td><td>Yes</td></tr> 
  <tr><td>filter_frequency</td><td></td></tr> 
  <tr><td>record_length</td><td>Not yet</td></tr> 
  <tr><td>sample_length</td> <td>Yes</td></tr>  
  <tr><td>gating</td><td>Not yet</td></tr>  
  <tr><td>gating_min</td><td>Not yet</td></tr>  
  <tr><td>gating_max</td><td>Not yet</td></tr>  
  <tr><td>probe_attenuation</td><td>Not yet</td></tr>  
  <tr><td>vertical_units</td><td>Not yet</td></tr>  
  <tr><td>vertical_offset</td><td>Not yet</td></tr>  
  <tr><td>vertical_scale</td><td>Not yet</td></tr>  
  <tr><td>vertical_position</td><td>Not yet</td></tr> 
  <tr><td>time</td><td>Yes</td></tr> 
  <tr><td>channels</td><td>Yes</td></tr>
<tr>
  <th colspan="2"> Optional </th>
</tr>
  <tr><td>Math</td><td>Yes</td></tr>
</table>
The channels variable is an array of channel objects that contains the name, vertical unit and the values of the used channels.

#### Functions
* fft(_scale('db')_,_window type_) 
* advancedFFT(_scale_,_window type_,_gatePosition_,_gateDuration_,_levelOffset_)
* plotChannel(_channel array_,_title array_,_xlimit array_) Creates
* plotMath(_channel array_,_title array_,_xlimit array_)

_For the list of all functions inside the scope file, please refer to the [header section in the scope.m](https://github.com/DimitriDeSchuyter/matlab-libraries/blob/master/scope/scope.m)_


### Eth class
PROFINET packets can be read
<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.










<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/DimitriDeSchuyter/matlab-libraries/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo.svg?style=for-the-badge
[forks-url]: https://github.com/DimitriDeSchuyter/matlab-libraries/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo.svg?style=for-the-badge
[stars-url]: https://github.com/DimitriDeSchuyter/matlab-libraries/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo.svg?style=for-the-badge
[issues-url]: https://github.com/DimitriDeSchuyter/matlab-libraries/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo.svg?style=for-the-badge
[license-url]: https://github.com/DimitriDeSchuyter/matlab-libraries/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/github_username

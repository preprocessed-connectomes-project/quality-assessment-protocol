template = 
'''
<!DOCTYPE html>
<html>
  <head>
    <style>
    body{
        margin: 0px;
    }
    ul {
        list-style-type: none;
        margin: 0;
        padding: 0;
        overflow: hidden;
        background-color: #333;
    }

    li {
        float: left;
    }

    li a {
        display: block;
        color: white;
        text-align: center;
        padding: 14px 16px;
        text-decoration: none;
    }

    li a:hover:not(.active) {
        background-color: #111;
    }

    .active {
        background-color: #4CAF50;
    }
    </style>
    <meta charset="UTF-8">
    <title>QAP REport for subject {subjectid}</title>
  </head>
  <body>
    <!-- start navbar -->
    <div class="navbar">
        <ul>
          <li><a href="#">{subjectid}</a></li>
          <li style="float:right"><a href="#about">QAP</a></li>
          <li style="float:right"><a href="#about">All Subjects</a></li>
          <li style="float:right"><a href="#about">Group Measures</a></li>
          
        </ul>
    </div>
    <!-- end navbar -->

    <!-- start Mean FD, DVARS, Global Signal -->
    <div id="meanfdplots">
        <h2>Mean FD, DVARS, Global Signal</h2>
    </div>
    <!-- end Mean FD, DVARS, Global Signal -->

    <!-- start Gray Plot -->
    <div id="grayplots">
        <h2>Gray Plot</h2>
    </div>
    <!-- end Gray Plot -->

    <!-- start Mean EPI Mosaic -->
    <div id="meanepi">
        <h2>Mean EPI Mosaic</h2>
    </div>
    <!-- end Mean EPI Mosaic -->

    <!-- start Signal Fluctuation Sensitivity Mosaic Mosaic -->
    <div id="sfs">
        <h2>Signal Fluctuation Sensitivity Mosaic</h2>
    </div>
    <!-- end Signal Fluctuation Sensitivity Mosaic Mosaic -->

  </body>
</html>
'''
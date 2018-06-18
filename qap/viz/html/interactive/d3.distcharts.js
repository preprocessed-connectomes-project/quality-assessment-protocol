/**
 * @fileOverview A D3 based distribution chart system. Supports: Box plots, Violin plots, Notched box plots, trend lines, beeswarm plot
 * @version 3.0
 */


/**
 * Creates a box plot, violin plot, and or notched box plot
 * @param settings Configuration options for the base plot
 * @param settings.data The data for the plot
 * @param settings.xName The name of the column that should be used for the x groups
 * @param settings.yName The name of the column used for the y values
 * @param {string} settings.selector The selector string for the main chart div
 * @param [settings.axisLabels={}] Defaults to the xName and yName
 * @param [settings.yTicks = 1] 1 = default ticks. 2 =  double, 0.5 = half
 * @param [settings.scale='linear'] 'linear' or 'log' - y scale of the chart
 * @param [settings.chartSize={width:800, height:400}] The height and width of the chart itself (doesn't include the container)
 * @param [settings.margin={top: 15, right: 60, bottom: 40, left: 50}] The margins around the chart (inside the main div)
 * @param [settings.constrainExtremes=false] Should the y scale include outliers?
 * @returns {object} chart A chart object
 */
function makeDistroChart(settings) {

    var chart = {};

    // Defaults
    chart.settings = {
        data: null,
        xName: null,
        yName: null,
        selector: null,
        axisLables: null,
        yTicks: 1,
        scale: 'linear',
        chartSize: {width: 800, height: 400},
        margin: {top: 15, right: 60, bottom: 40, left: 50},
        constrainExtremes: false,
        color: d3.scale.category10()
    };
    for (var setting in settings) {
        chart.settings[setting] = settings[setting]
    }


    function formatAsFloat(d) {
        if (d % 1 !== 0) {
            return d3.format(".2f")(d);
        } else {
            return d3.format(".0f")(d);
        }
    }

    function logFormatNumber(d) {
        var x = Math.log(d) / Math.log(10) + 1e-6;
        return Math.abs(x - Math.floor(x)) < 0.6 ? formatAsFloat(d) : "";
    }

    chart.yFormatter = formatAsFloat;

    chart.data = chart.settings.data;

    chart.groupObjs = {}; //The data organized by grouping and sorted as well as any metadata for the groups
    chart.objs = {mainDiv: null, chartDiv: null, g: null, xAxis: null, yAxis: null};
    chart.colorFunct = null;

    /**
     * Takes an array, function, or object mapping and created a color function from it
     * @param {function|[]|object} colorOptions
     * @returns {function} Function to be used to determine chart colors
     */
    function getColorFunct(colorOptions) {
        if (typeof colorOptions == 'function') {
            return colorOptions
        } else if (Array.isArray(colorOptions)) {
            //  If an array is provided, map it to the domain
            var colorMap = {}, cColor = 0;
            for (var cName in chart.groupObjs) {
                colorMap[cName] = colorOptions[cColor];
                cColor = (cColor + 1) % colorOptions.length;
            }
            return function (group) {
                return colorMap[group];
            }
        } else if (typeof colorOptions == 'object') {
            // if an object is provided, assume it maps to  the colors
            return function (group) {
                return colorOptions[group];
            }
        } else {
            return d3.scale.category10();
        }
    }

    /**
     * Takes a percentage as returns the values that correspond to that percentage of the group range witdh
     * @param objWidth Percentage of range band
     * @param gName The bin name to use to get the x shift
     * @returns {{left: null, right: null, middle: null}}
     */
    function getObjWidth(objWidth, gName) {
        var objSize = {left: null, right: null, middle: null};
        var width = chart.xScale.rangeBand() * (objWidth / 100);
        var padding = (chart.xScale.rangeBand() - width) / 2;
        var gShift = chart.xScale(gName);
        objSize.middle = chart.xScale.rangeBand() / 2 + gShift;
        objSize.left = padding + gShift;
        objSize.right = objSize.left + width;
        return objSize;
    }

    /**
     * Adds jitter to the  scatter point plot
     * @param doJitter true or false, add jitter to the point
     * @param width percent of the range band to cover with the jitter
     * @returns {number}
     */
    function addJitter(doJitter, width) {
        if (doJitter !== true || width == 0) {
            return 0
        }
        return Math.floor(Math.random() * width) - width / 2;
    }

    function shallowCopy(oldObj) {
        var newObj = {};
        for (var i in oldObj) {
            if (oldObj.hasOwnProperty(i)) {
                newObj[i] = oldObj[i];
            }
        }
        return newObj;
    }

    /**
     * Closure that creates the tooltip hover function
     * @param groupName Name of the x group
     * @param metrics Object to use to get values for the group
     * @returns {Function} A function that provides the values for the tooltip
     */
    function tooltipHover(groupName, group) {
        var tooltipString = "";
        tooltipString += "Max: " + formatAsFloat(group.metrics.max, 0.1);
        tooltipString += "<br />Q3: " + formatAsFloat(group.metrics.quartile3);
        tooltipString += "<br />Median: " + formatAsFloat(group.metrics.median);
        tooltipString += "<br />Q1: " + formatAsFloat(group.metrics.quartile1);
        tooltipString += "<br />Min: " + formatAsFloat(group.metrics.min);

        if (group.marks && group.marks.length > 0) {
            tooltipString += "<hr />Marks:<ul>";
            for (var i in group.marks) {
                tooltipString += "<li>" + group.marks[i].toFixed(2) + "</li>"
            }
            tooltipString += "</ul>"
        }

        return function () {
            chart.objs.tooltip.transition().duration(200).style("opacity", 1.0);
            chart.objs.tooltip.html(tooltipString)
        };
    }

    /**
     * Parse the data and calculates base values for the plots
     */
    !function prepareData() {
        function calcMetrics(values) {

            var metrics = { //These are the original non�scaled values
                max: null,
                upperOuterFence: null,
                upperInnerFence: null,
                quartile3: null,
                median: null,
                mean: null,
                iqr: null,
                quartile1: null,
                lowerInnerFence: null,
                lowerOuterFence: null,
                min: null
            };

            metrics.min = d3.min(values);
            metrics.quartile1 = d3.quantile(values, 0.25);
            metrics.median = d3.median(values);
            metrics.mean = d3.mean(values);
            metrics.quartile3 = d3.quantile(values, 0.75);
            metrics.max = d3.max(values);
            metrics.iqr = metrics.quartile3 - metrics.quartile1;

            //The inner fences are the closest value to the IQR without going past it (assumes sorted lists)
            var LIF = metrics.quartile1 - (1.5 * metrics.iqr);
            var UIF = metrics.quartile3 + (1.5 * metrics.iqr);
            for (var i = 0; i <= values.length; i++) {
                if (values[i] < LIF) {
                    continue;
                }
                if (!metrics.lowerInnerFence && values[i] >= LIF) {
                    metrics.lowerInnerFence = values[i];
                    continue;
                }
                if (values[i] > UIF) {
                    metrics.upperInnerFence = values[i - 1];
                    break;
                }
            }

            metrics.lowerOuterFence = metrics.quartile1 - (3 * metrics.iqr);
            metrics.upperOuterFence = metrics.quartile3 + (3 * metrics.iqr);
            if (!metrics.lowerInnerFence) {
                metrics.lowerInnerFence = metrics.min;
            }
            if (!metrics.upperInnerFence) {
                metrics.upperInnerFence = metrics.max;
            }
            return metrics
        }

        var current_x = null;
        var current_y = null;
        var current_row;

        // Group the values
        for (current_row = 0; current_row < chart.data.length; current_row++) {
            current_x = chart.data[current_row][chart.settings.xName];
            current_y = chart.data[current_row][chart.settings.yName];
            current_isMark = !!chart.data[current_row].mark

            if (chart.groupObjs.hasOwnProperty(current_x)) {
                chart.groupObjs[current_x].values.push(current_y);
            } else {
                chart.groupObjs[current_x] = { marks: [], values: [current_y] };
            }
            if (current_isMark) {
                chart.groupObjs[current_x].marks.push(current_y)
            }
        }

        for (var cName in chart.groupObjs) {
            chart.groupObjs[cName].values.sort(d3.ascending);
            chart.groupObjs[cName].marks.sort(d3.ascending);
            chart.groupObjs[cName].metrics = {};
            chart.groupObjs[cName].metrics = calcMetrics(chart.groupObjs[cName].values);

        }
    }();

    /**
     * Prepare the chart settings and chart div and svg
     */
    !function prepareSettings() {
        //Set base settings
        chart.margin = chart.settings.margin;
        chart.divWidth = chart.settings.chartSize.width;
        chart.divHeight = chart.settings.chartSize.height;
        chart.width = chart.divWidth - chart.margin.left - chart.margin.right;
        chart.height = chart.divHeight - chart.margin.top - chart.margin.bottom;

        if (chart.settings.axisLabels) {
            chart.xAxisLable = chart.settings.axisLabels.xAxis;
            chart.yAxisLable = chart.settings.axisLabels.yAxis;
        } else {
            chart.xAxisLable = chart.settings.xName;
            chart.yAxisLable = chart.settings.yName;
        }

        if (chart.settings.scale === 'log') {
            chart.yScale = d3.scale.log();
            chart.yFormatter = logFormatNumber;
        } else {
            chart.yScale = d3.scale.linear();
        }

        if (chart.settings.constrainExtremes === true) {
            var fences = [];
            for (var cName in chart.groupObjs) {
                fences.push(chart.groupObjs[cName].metrics.lowerInnerFence);
                fences.push(chart.groupObjs[cName].metrics.upperInnerFence);
            }
            chart.range = d3.extent(fences);

        } else {
            chart.range = d3.extent(chart.data, function (d) {return d[chart.settings.yName];});
        }

        chart.colorFunct = getColorFunct(chart.settings.colors);

        // Build Scale functions
        chart.yScale.range([chart.height, 0]).domain(chart.range).nice().clamp(true);
        chart.xScale = d3.scale.ordinal().domain(Object.keys(chart.groupObjs)).rangeBands([0, chart.width]);

        //Build Axes Functions
        chart.objs.yAxis = d3.svg.axis()
            .scale(chart.yScale)
            .orient("left")
            .tickFormat(chart.yFormatter)
            .outerTickSize(0)
            .innerTickSize(-chart.width + (chart.margin.right + chart.margin.left));
        chart.objs.yAxis.ticks(chart.objs.yAxis.ticks()*chart.settings.yTicks);
        chart.objs.xAxis = d3.svg.axis().scale(chart.xScale).orient("bottom").tickSize(5);
    }();

    /**
     * Updates the chart based on the current settings and window size
     * @returns {*}
     */
    chart.update = function () {
        // Update chart size based on view port size
        chart.width = parseInt(chart.objs.chartDiv.style("width"), 10) - (chart.margin.left + chart.margin.right);
        chart.height = parseInt(chart.objs.chartDiv.style("height"), 10) - (chart.margin.top + chart.margin.bottom);

        // Update scale functions
        chart.xScale.rangeBands([0, chart.width]);
        chart.yScale.range([chart.height, 0]);

        // Update the yDomain if the Violin plot clamp is set to -1 meaning it will extend the violins to make nice points
        if (chart.violinPlots && chart.violinPlots.options.show == true && chart.violinPlots.options._yDomainVP != null) {
            chart.yScale.domain(chart.violinPlots.options._yDomainVP).nice().clamp(true);
        } else {
            chart.yScale.domain(chart.range).nice().clamp(true);
        }

        //Update axes
        chart.objs.g.select('.x.axis').attr("transform", "translate(0," + chart.height + ")").call(chart.objs.xAxis)
            .selectAll("text")
            .attr("y", 20)
            .attr("x", -5)
            .style("text-anchor", "middle");
        chart.objs.g.select('.x.axis .label').attr("x", chart.width / 2);
        chart.objs.g.select('.y.axis').call(chart.objs.yAxis.innerTickSize(-chart.width));
        chart.objs.g.select('.y.axis .label').attr("x", -chart.height / 2);
        chart.objs.chartDiv.select('svg').attr("width", chart.width + (chart.margin.left + chart.margin.right)).attr("height", chart.height + (chart.margin.top + chart.margin.bottom));

        return chart;
    };

    /**
     * Prepare the chart html elements
     */
    !function prepareChart() {
        // Build main div and chart div
        chart.objs.mainDiv = d3.select(chart.settings.selector)
            .style("max-width", chart.divWidth + "px");
        // Add all the divs to make it centered and responsive
        chart.objs.mainDiv.append("div")
            .attr("class", "inner-wrapper")
            .style("padding-bottom", (chart.divHeight / chart.divWidth) * 100 + "%")
            .append("div").attr("class", "outer-box")
            .append("div").attr("class", "inner-box");
        // Capture the inner div for the chart (where the chart actually is)
        chart.selector = chart.settings.selector + " .inner-box";
        chart.objs.chartDiv = d3.select(chart.selector);
        d3.select(window).on('resize.' + chart.selector, chart.update);

        // Create the svg
        chart.objs.g = chart.objs.chartDiv.append("svg")
            .attr("class", "chart-area")
            .attr("width", chart.width + (chart.margin.left + chart.margin.right))
            .attr("height", chart.height + (chart.margin.top + chart.margin.bottom))
            .append("g")
            .attr("transform", "translate(" + chart.margin.left + "," + chart.margin.top + ")");

        // Create axes
        chart.objs.axes = chart.objs.g.append("g").attr("class", "axis");
        chart.objs.axes.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + chart.height + ")")
            .call(chart.objs.xAxis);
        chart.objs.axes.append("g")
            .attr("class", "y axis")
            .call(chart.objs.yAxis)
            .append("text")
            .attr("class", "label")
            .attr("transform", "rotate(-90)")
            .attr("y", -42)
            .attr("x", -chart.height / 2)
            .attr("dy", ".71em")
            .style("text-anchor", "middle")
            .text(chart.yAxisLable);

        // Create tooltip div
        chart.objs.tooltip = chart.objs.mainDiv.append('div').attr('class', 'tooltip');
        for (var cName in chart.groupObjs) {
            chart.groupObjs[cName].g = chart.objs.g.append("g").attr("class", "group");
            chart.groupObjs[cName].g.on("mouseover", function () {
                chart.objs.tooltip
                    .style("display", null)
                    .style("left", (d3.event.pageX) + "px")
                    .style("top", (d3.event.pageY - 28) + "px");
            }).on("mouseout", function () {
                chart.objs.tooltip.style("display", "none");
            }).on("mousemove", tooltipHover(cName, chart.groupObjs[cName]))
        }
        chart.update();
    }();

    /**
     * Render a box plot on the current chart
     * @param options
     * @param [options.show=true] Toggle the whole plot on and off
     * @param [options.showBox=true] Show the box part of the box plot
     * @param [options.showWhiskers=true] Show the whiskers
     * @param [options.showMedian=true] Show the median line
     * @param [options.showMean=false] Show the mean line
     * @param [options.medianCSize=3] The size of the circle on the median
     * @param [options.showOutliers=true] Plot outliers
     * @param [options.boxwidth=30] The max percent of the group rangeBand that the box can be
     * @param [options.lineWidth=boxWidth] The max percent of the group rangeBand that the line can be
     * @param [options.outlierScatter=false] Spread out the outliers so they don't all overlap (in development)
     * @param [options.outlierCSize=2] Size of the outliers
     * @param [options.colors=chart default] The color mapping for the box plot
     * @returns {*} The chart object
     */
    chart.renderBoxPlot = function (options) {
        chart.boxPlots = {};

        // Defaults
        var defaultOptions = {
            show: true,
            showBox: true,
            showWhiskers: true,
            showMedian: true,
            showMean: false,
            medianCSize: 3.5,
            showOutliers: true,
            boxWidth: 30,
            lineWidth: null,
            scatterOutliers: false,
            outlierCSize: 2.5,
            colors: chart.colorFunct
        };
        chart.boxPlots.options = shallowCopy(defaultOptions);
        for (var option in options) {
            chart.boxPlots.options[option] = options[option]
        }
        var bOpts = chart.boxPlots.options;

        //Create box plot objects
        for (var cName in chart.groupObjs) {
            chart.groupObjs[cName].boxPlot = {};
            chart.groupObjs[cName].boxPlot.objs = {};
        }


        /**
         * Calculates all the outlier points for each group
         */
        !function calcAllOutliers() {

            /**
             * Create lists of the outliers for each content group
             * @param cGroup The object to modify
             * @return null Modifies the object in place
             */
            function calcOutliers(cGroup) {
                var cExtremes = [];
                var cOutliers = [];
                var cOut, idx;
                for (idx = 0; idx <= cGroup.values.length; idx++) {
                    cOut = {value: cGroup.values[idx]};
                    if (cOut.value < cGroup.metrics.lowerInnerFence) {
                        if (cOut.value < cGroup.metrics.lowerOuterFence) {
                            cExtremes.push(cOut);
                        } else {
                            cOutliers.push(cOut);
                        }
                    } else if (cOut.value > cGroup.metrics.upperInnerFence) {
                        if (cOut.value > cGroup.metrics.upperOuterFence) {
                            cExtremes.push(cOut);
                        } else {
                            cOutliers.push(cOut);
                        }
                    }
                }
                cGroup.boxPlot.objs.outliers = cOutliers;
                cGroup.boxPlot.objs.extremes = cExtremes;

                cGroup.boxPlot.objs.marks = []
                for (idx = 0; idx < cGroup.marks.length; idx++) {
                    cOut = {value: cGroup.marks[idx]};
                    cGroup.boxPlot.objs.marks.push(cOut)
                }
            }

            for (var cName in chart.groupObjs) {
                calcOutliers(chart.groupObjs[cName]);
            }
        }();

        /**
         * Take updated options and redraw the box plot
         * @param updateOptions
         */
        chart.boxPlots.change = function (updateOptions) {
            if (updateOptions) {
                for (var key in updateOptions) {
                    bOpts[key] = updateOptions[key]
                }
            }

            for (var cName in chart.groupObjs) {
                chart.groupObjs[cName].boxPlot.objs.g.remove()
            }
            chart.boxPlots.prepareBoxPlot();
            chart.boxPlots.update()
        };

        chart.boxPlots.reset = function () {
            chart.boxPlots.change(defaultOptions)
        };
        chart.boxPlots.show = function (opts) {
            if (opts !== undefined) {
                opts.show = true;
                if (opts.reset) {
                    chart.boxPlots.reset()
                }
            } else {
                opts = {show: true};
            }
            chart.boxPlots.change(opts)

        };
        chart.boxPlots.hide = function (opts) {
            if (opts !== undefined) {
                opts.show = false;
                if (opts.reset) {
                    chart.boxPlots.reset()
                }
            } else {
                opts = {show: false};
            }
            chart.boxPlots.change(opts)
        };

        /**
         * Update the box plot obj values
         */
        chart.boxPlots.update = function () {
            var cName, cBoxPlot;

            for (cName in chart.groupObjs) {
                cBoxPlot = chart.groupObjs[cName].boxPlot;

                // Get the box width
                var objBounds = getObjWidth(bOpts.boxWidth, cName);
                var width = (objBounds.right - objBounds.left);

                var sMetrics = {}; //temp var for scaled (plottable) metric values
                for (var attr in chart.groupObjs[cName].metrics) {
                    sMetrics[attr] = null;
                    sMetrics[attr] = chart.yScale(chart.groupObjs[cName].metrics[attr]);
                }

                // Box
                if (cBoxPlot.objs.box) {
                    cBoxPlot.objs.box
                        .attr("x", objBounds.left)
                        .attr('width', width)
                        .attr("y", sMetrics.quartile3)
                        .attr("rx", 1)
                        .attr("ry", 1)
                        .attr("height", -sMetrics.quartile3 + sMetrics.quartile1)
                }

                // Lines
                var lineBounds = null;
                if (bOpts.lineWidth) {
                    lineBounds = getObjWidth(bOpts.lineWidth, cName)
                } else {
                    lineBounds = objBounds
                }
                // --Whiskers
                if (cBoxPlot.objs.upperWhisker) {
                    cBoxPlot.objs.upperWhisker.fence
                        .attr("x1", lineBounds.left)
                        .attr("x2", lineBounds.right)
                        .attr('y1', sMetrics.upperInnerFence)
                        .attr("y2", sMetrics.upperInnerFence);
                    cBoxPlot.objs.upperWhisker.line
                        .attr("x1", lineBounds.middle)
                        .attr("x2", lineBounds.middle)
                        .attr('y1', sMetrics.quartile3)
                        .attr("y2", sMetrics.upperInnerFence);

                    cBoxPlot.objs.lowerWhisker.fence
                        .attr("x1", lineBounds.left)
                        .attr("x2", lineBounds.right)
                        .attr('y1', sMetrics.lowerInnerFence)
                        .attr("y2", sMetrics.lowerInnerFence);
                    cBoxPlot.objs.lowerWhisker.line
                        .attr("x1", lineBounds.middle)
                        .attr("x2", lineBounds.middle)
                        .attr('y1', sMetrics.quartile1)
                        .attr("y2", sMetrics.lowerInnerFence);
                }

                // --Median
                if (cBoxPlot.objs.median) {
                    cBoxPlot.objs.median.line
                        .attr("x1", lineBounds.left)
                        .attr("x2", lineBounds.right)
                        .attr('y1', sMetrics.median)
                        .attr("y2", sMetrics.median);
                    cBoxPlot.objs.median.circle
                        .attr("cx", lineBounds.middle)
                        .attr("cy", sMetrics.median)
                }

                // --Mean
                if (cBoxPlot.objs.mean) {
                    cBoxPlot.objs.mean.line
                        .attr("x1", lineBounds.left)
                        .attr("x2", lineBounds.right)
                        .attr('y1', sMetrics.mean)
                        .attr("y2", sMetrics.mean);
                    cBoxPlot.objs.mean.circle
                        .attr("cx", lineBounds.middle)
                        .attr("cy", sMetrics.mean);
                }

                // Outliers
                var pt;
                if (cBoxPlot.objs.outliers) {
                    for (pt in cBoxPlot.objs.outliers) {
                        cBoxPlot.objs.outliers[pt].point
                            .attr("cx", objBounds.middle + addJitter(bOpts.scatterOutliers, width))
                            .attr("cy", chart.yScale(cBoxPlot.objs.outliers[pt].value));
                    }
                }
                if (cBoxPlot.objs.extremes) {
                    for (pt in cBoxPlot.objs.extremes) {
                        cBoxPlot.objs.extremes[pt].point
                            .attr("cx", objBounds.middle + addJitter(bOpts.scatterOutliers, width))
                            .attr("cy", chart.yScale(cBoxPlot.objs.extremes[pt].value));
                    }
                }
                if (cBoxPlot.objs.marks) {
                    for (pt in cBoxPlot.objs.marks) {
                        cBoxPlot.objs.marks[pt].point
                            .attr("cx", objBounds.middle + addJitter(bOpts.scatterOutliers, width))
                            .attr("cy", chart.yScale(cBoxPlot.objs.marks[pt].value));
                    }
                }
            }
        };

        /**
         * Create the svg elements for the box plot
         */
        chart.boxPlots.prepareBoxPlot = function () {
            var cName, cBoxPlot;

            if (bOpts.colors) {
                chart.boxPlots.colorFunct = getColorFunct(bOpts.colors);
            } else {
                chart.boxPlots.colorFunct = chart.colorFunct
            }

            if (bOpts.show == false) {
                return
            }

            for (cName in chart.groupObjs) {
                cBoxPlot = chart.groupObjs[cName].boxPlot;

                cBoxPlot.objs.g = chart.groupObjs[cName].g.append("g").attr("class", "box-plot");

                //Plot Box (default show)
                if (bOpts.showBox) {
                    cBoxPlot.objs.box = cBoxPlot.objs.g.append("rect")
                        .attr("class", "box")
                        .style("fill", chart.boxPlots.colorFunct(cName))
                        .style("stroke", chart.boxPlots.colorFunct(cName));
                    //A stroke is added to the box with the group color, it is
                    // hidden by default and can be shown through css with stroke-width
                }

                //Plot Median (default show)
                if (bOpts.showMedian) {
                    cBoxPlot.objs.median = {line: null, circle: null};
                    cBoxPlot.objs.median.line = cBoxPlot.objs.g.append("line")
                        .attr("class", "median");
                    cBoxPlot.objs.median.circle = cBoxPlot.objs.g.append("circle")
                        .attr("class", "median")
                        .attr('r', bOpts.medianCSize)
                        .style("fill", chart.boxPlots.colorFunct(cName));
                }

                // Plot Mean (default no plot)
                if (bOpts.showMean) {
                    cBoxPlot.objs.mean = {line: null, circle: null};
                    cBoxPlot.objs.mean.line = cBoxPlot.objs.g.append("line")
                        .attr("class", "mean");
                    cBoxPlot.objs.mean.circle = cBoxPlot.objs.g.append("circle")
                        .attr("class", "mean")
                        .attr('r', bOpts.medianCSize)
                        .style("fill", chart.boxPlots.colorFunct(cName));
                }

                // Plot Whiskers (default show)
                if (bOpts.showWhiskers) {
                    cBoxPlot.objs.upperWhisker = {fence: null, line: null};
                    cBoxPlot.objs.lowerWhisker = {fence: null, line: null};
                    cBoxPlot.objs.upperWhisker.fence = cBoxPlot.objs.g.append("line")
                        .attr("class", "upper whisker")
                        .style("stroke", chart.boxPlots.colorFunct(cName));
                    cBoxPlot.objs.upperWhisker.line = cBoxPlot.objs.g.append("line")
                        .attr("class", "upper whisker")
                        .style("stroke", chart.boxPlots.colorFunct(cName));

                    cBoxPlot.objs.lowerWhisker.fence = cBoxPlot.objs.g.append("line")
                        .attr("class", "lower whisker")
                        .style("stroke", chart.boxPlots.colorFunct(cName));
                    cBoxPlot.objs.lowerWhisker.line = cBoxPlot.objs.g.append("line")
                        .attr("class", "lower whisker")
                        .style("stroke", chart.boxPlots.colorFunct(cName));
                }

                // Plot outliers (default show)
                if (bOpts.showOutliers) {
                    if (!cBoxPlot.objs.outliers) calcAllOutliers();
                    var pt;
                    if (cBoxPlot.objs.outliers.length) {
                        var outDiv = cBoxPlot.objs.g.append("g").attr("class", "boxplot outliers");
                        for (pt in cBoxPlot.objs.outliers) {
                            cBoxPlot.objs.outliers[pt].point = outDiv.append("circle")
                                .attr("class", "outlier")
                                .attr('r', bOpts.outlierCSize)
                                .style("fill", chart.boxPlots.colorFunct(cName));
                        }
                    }

                    if (cBoxPlot.objs.extremes.length) {
                        var extDiv = cBoxPlot.objs.g.append("g").attr("class", "boxplot extremes");
                        for (pt in cBoxPlot.objs.extremes) {
                            cBoxPlot.objs.extremes[pt].point = extDiv.append("circle")
                                .attr("class", "extreme")
                                .attr('r', bOpts.outlierCSize)
                                .style("stroke", chart.boxPlots.colorFunct(cName));
                        }
                    }

                    if (cBoxPlot.objs.marks.length) {
                        var extDiv = cBoxPlot.objs.g.append("g").attr("class", "boxplot marks");
                        for (pt in cBoxPlot.objs.marks) {
                            cBoxPlot.objs.marks[pt].point = extDiv.append("circle")
                                .attr("class", "mark")
                                .attr('r', bOpts.outlierCSize)
                                .style("stroke", chart.boxPlots.colorFunct(cName));
                        }
                    }
                }


            }
        };
        chart.boxPlots.prepareBoxPlot();

        d3.select(window).on('resize.' + chart.selector + '.boxPlot', chart.boxPlots.update);
        chart.boxPlots.update();
        return chart;

    };

    return chart;
}
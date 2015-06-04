var phylocanvas = undefined;

HTMLWidgets.widget({

  name: 'radial_phylo',

  type: 'output',

  initialize: function(el, width, height) {

    return {
        width: el.offsetWidth,
        height: el.offsetHeight
    }

  },

  renderValue: function(el, x, instance) {

    // function to add css to <head>
    var addCSS = function(css) {
        var head = document.head || document.getElementsByTagName('head')[0];
        var style = document.createElement('style');
        style.type = "text/css";
        if (style.styleSheet) {
            style.styleSheet.cssText = css;
        } else {
            style.appendChild(document.createTextNode(css));
        }
        head.appendChild(style);
    }
                   
    // if 'canvas_size' is a number then use it
    var width = undefined;
    var height = undefined;
    if (typeof(x.canvas_size === 'number') && x.canvas_size%1 === 0) {
        width = x.canvas_size;
        height = x.canvas_size;
    } else {
        // get width and height of current window
        width = Math.max(instance.width, instance.height);
        height = Math.max(instance.width, instance.height);
    }

    // if user provides large canvas_size give the svg that width and height
    if (x.scale === true) {
        addCSS("svg { width: 100%; height: 100%; }");
    } else {
        addCSS("svg { height: " + width + "px; width: " + height + "px; }");
        addCSS("body { overflow: scroll !important; }");
    }

    // set some options of the phylogram
    Smits.PhyloCanvas.Render.Parameters.Circular['bufferRadius'] = .42;
    Smits.PhyloCanvas.Render.Parameters.Circular['bufferOuterLabels'] = 0;

    // get xml file's contents and plot the phylogram
    var xml_relpath = HTMLWidgets.getAttachmentUrl('phyloxml', 'xml');
    $.get(xml_relpath, function(xmldata) {

        var dataObject = {
            phyloxml: xmldata,
            fileSource: true
        };

        phylocanvas = new Smits.PhyloCanvas(
            dataObject,
            el.id,
            width, height,
            'circular',
            x.scale  // holds true or false
        );

    });

	// Create 'Save Image' text link
	var a = document.createElement('a');
	var downloadLink = document.createTextNode("Save Image (Web Browser Only)");
	a.appendChild(downloadLink);
	a.href = "#"
	a.title = "Save Image";
    a.id = "download_link";
	var widget = document.body.children[0];
	document.body.insertBefore(a, widget);

	// Attach click handler to save the image when link clicked
	a.onclick = function() {
		$(document.body).append("<canvas id='canvg'></canvas>");
		$(document.body).append("<img id='svgimg' src=''>");
        // The anchor below allows you to suggest a filename for the image. You
        // have to download the image through it, by setting the href to the
        // image src once it's ready and setting the filename in the 'downloaad'
        // attribute.
        var download_anchor = document.createElement("a");
		
		// Use Raphael.Export to get the SVG from the phylogram
		var s = phylocanvas.getSvg();
		var c = s.svg;
		var svg = c.toSVG();
		
		// Use canvg to draw the SVG onto the empty canvas
		canvg(document.getElementById("canvg"), svg);
		
		// Give the canvas a second or two to load, then set the image source
		setTimeout(function() {
		           //fetch the dataURL from the canvas and set it as src on the image
		           var dataURL = document.getElementById("canvg").toDataURL("image/png");
		           document.getElementById("svgimg").src = dataURL;
        }, 1000);

        // Initiate the download
        setTimeout(function() {
            var svgimg = document.getElementById("svgimg");
            // Indicate the data is a byte stream so it doesn't open image in browser
            var img_data_uri =  svgimg.src.replace("image/png", "image/octet-stream");
            $(download_anchor).attr("download", "phylo.png");
            download_anchor.href = img_data_uri;
            // Fake a mouse click on our anchor to initiate the download
            if (document.createEvent) {
                var e = document.createEvent("MouseEvents");
                e.initMouseEvent("click", true, true, window,
                    0, 0, 0, 0, 0, false, false, false,
                    false, 0, null);
                download_anchor.dispatchEvent(e);
            } else if (download_anchor.fireEvent) {
                download_anchor.fireEvent("onclick");
            }
        }, 1000);
    }

  },

  resize: function(el, width, height, instance) {

    // maximize the size of the phylogram by
    var size = Math.max(width, height);

    // delete the existing phylogram and remove its svg from the dom
    ////delete phylocanvas;
    ////var svg = el.firstElementChild;
    ////el.removeChild(svg);

    ////
    ////el.innerText = x;

    // phylocanvas = new Smits.PhyloCanvas(
    //     {
    //         newick: x.dataObject
    //     },
    //     el.id,
    //     size, size,
    //     'circular'
    // );

  }

});

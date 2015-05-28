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
                   
    // if 'radius' is a number then use it
    var width = undefined;
    var height = undefined;
    if (typeof(x.radius === 'number') && x.radius%1 === 0) {
        width = x.radius;
        height = x.radius;
    } else {
        // get width and height of current window
        width = Math.max(instance.width, instance.height);
        height = Math.max(instance.width, instance.height);
    }

    // if user provides large radius give the svg that width and height
    if (x.autoResize === true) {
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

        var phylocanvas = new Smits.PhyloCanvas(
            dataObject,
            el.id,
            width, height,
            'circular',
            x.autoResize  // holds true or false
        );

    });

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

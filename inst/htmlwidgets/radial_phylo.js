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

    // add css
    var css = "svg { height: 100%; width: 100%; }";
    var head = document.head || document.getElementsByTagName('head')[0];
    var style = document.createElement('style');
    style.type = "text/css";
    if (style.styleSheet) {
        style.styleSheet.cssText = css;
    } else {
        style.appendChild(document.createTextNode(css));
    }
    head.appendChild(style);
                   
    // if 'radius' is a number then use it
    var width = undefined;
    var height = undefined;
    if (typeof(x.radius === 'number') && x.radius%1 === 0) {
        width = x.radius;
        height = x.radius;
    } else {
        // get width and height of current window
        width = el.offsetWidth;
        height = el.offsetHeight;
    }

    // set some options of the phylogram
    Smits.PhyloCanvas.Render.Parameters.Circular['bufferRadius'] = .35;
    Smits.PhyloCanvas.Render.Parameters.Circular['bufferOuterLabels'] = 0;

    // plot the phylgram
    var phylocanvas = new Smits.PhyloCanvas(
        {
            newick: x.dataObject
        },
        el.id,
        width, height,
        'circular',
        x.autoResize  // holds true or false
    );

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

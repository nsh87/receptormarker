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

    // get width and height of current window
    var width = el.offsetWidth;
    var height = el.offsetHeight;

    // set some options of the phylogram
    Smits.PhyloCanvas.Render.Parameters.Circular['bufferRadius'] = .35;

    // plot the phylgram
    var phylocanvas = new Smits.PhyloCanvas(
        {
            newick: x.dataObject
        },
        el.id,
        width, height,
        'circular'
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

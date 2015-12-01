HTMLWidgets.widget({

  name: 'convergence',

  type: 'output',

  initialize: function(el, width, height) {

    return {
      width: el.offsetWidth,
      height: el.offsetHeight
    }

  },

  renderValue: function(el, x, instance) {
    var xml_path = HTMLWidgets.getAttachmentUrl('convergencexml','file');
    
    var request = new XMLHttpRequest();
    request.open("GET", xml_path, false);
    request.send();
            console.log(request);
    var network_xml = request.responseText;
    
    var network = {
        // you need to specify a data schema for custom attributes!
        dataSchema: {
            nodes: [ { name: "label", type: "string" },
                { name: "foo", type: "string" }
            ],
            edges: [ { name: "label", type: "string" },
                { name: "bar", type: "string" }
            ]
        },
        // NOTE the custom attributes on nodes and edges
        data: {
            nodes: [ { id: "1", label: "1", foo: "Is this the real life?" },
                { id: "2", label: "2", foo: "Is this just fantasy?" }
            ],
            edges: [ { id: "2to1", target: "1", source: "2", label: "2 to 1", bar: "Caught in a landslide..." }
            ]
        }
    };
    
    var options = {
      swfPath: HTMLWidgets.getAttachmentUrl('cytoscapeweb', 'CytoscapeWeb'),
      flashInstallerPath: HTMLWidgets.getAttachmentUrl('cytoscapeweb', 'playerProductInstall')
    };
    var vis_style = {
      global: {
          backgroundColor: "#ABCFD6"
      },
      nodes: {
          shape: "CIRCLE",
          borderWidth: 3,
          borderColor: "#ffffff",
          color: "#0b94b1",
          size: 25,
          labelHorizontalAnchor: "center"
      },
      edges: {
          width: 3,
          color: "#0B94B1"
      }
    };
    var draw_options = {

      //network: network,
      network: network_xml,

      nodeLabelsVisible: x.isLabel,

      // let's try another layout
      layout: "Circle",

      // set the style at initialisation
      visualStyle: vis_style,

      // hide pan zoom
      panZoomControlVisible: false
    };
    console.log(JSON.stringify(draw_options));
    instance.cy = new org.cytoscapeweb.Visualization(el.id, options);
    instance.cy.draw(draw_options);
  },

  resize: function(el, width, height, instance) {

  }

});

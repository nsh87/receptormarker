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
    
    window.onload = function() {
      
      var options = {
        swfPath: HTMLWidgets.getAttachmentUrl('cytoscapeweb', 'CytoscapeWeb'),
        flashInstallerPath: HTMLWidgets.getAttachmentUrl('cytoscapeweb',
                                                         'playerProductInstall')
      };
      
      var vis_style = {
        global: {
            backgroundColor: x.background_color
        },
        nodes: {
            shape: x.node_shape,
            borderWidth: x.border_width,
            borderColor: x.border_color,
            color: x.node_color,
            size: x.node_size,
            labelHorizontalAnchor: x.label_horizontal_pos,
            labelVerticalAnchor: x.label_vertical_pos
        },
        edges: {
            width: x.edge_width,
            color: x.edge_color,
            targetArrowShape: "NONE"
        }
      };
      
      // Draw options at http://cytoscapeweb.cytoscape.org/documentation
      var draw_options = {
  
        //network: network,
        network: x.xml_string,
        
        //nodeLabelsVisible: x.isLabel,
        nodeLabelsVisible: true,
  
        /* A layout object or name: http://cytoscapeweb.cytoscape.org/
           documentation/layout#section/Layout */
        layout: {
          name: "ForceDirected",
          options: {
            weightAttr: "weight",
            autoStabilize: true,
            maxTime: 2000
          }
        },
  
        // Set the style at initialisation
        visualStyle: vis_style,
        
        // Show pan zoom
        panZoomControlVisible: false,
        
        panZoomControlPosition: "bottomRight"
        
      };
      instance.cy = new org.cytoscapeweb.Visualization(el.id, options);
      instance.cy.draw(draw_options);
      
    };
    
  },

  resize: function(el, width, height, instance) {
    
    instance.cy.zoomToFit();

  }

});

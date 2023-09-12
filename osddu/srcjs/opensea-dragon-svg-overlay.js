// OpenSeadragon SVG Overlay plugin 0.0.5

(function() {

    var $ = window.OpenSeadragon;

    if (!$) {
        $ = require('openseadragon');
        if (!$) {
            throw new Error('OpenSeadragon is missing.');
        }
    }

    var svgNS = 'http://www.w3.org/2000/svg';

    // ----------
    $.Viewer.prototype.svgOverlay = function() {
        if (this._svgOverlayInfo) {
            return this._svgOverlayInfo;
        }

        this._svgOverlayInfo = new Overlay(this);
        return this._svgOverlayInfo;
    };

    // ----------
    var Overlay = function(viewer) {
        var self = this;

        this._viewer = viewer;
        this._containerWidth = 0;
        this._containerHeight = 0;

        this._svg = document.createElementNS(svgNS, 'svg');
        this._svg.style.position = 'absolute';
        this._svg.style.left = 0;
        this._svg.style.top = 0;
        this._svg.style.width = '100%';
        this._svg.style.height = '100%';
        this._viewer.canvas.appendChild(this._svg);

        this._node = document.createElementNS(svgNS, 'g');
        this._svg.appendChild(this._node);

        this._viewer.addHandler('animation', function() {
            self.resize();
        });

        this._viewer.addHandler('open', function() {
            self.resize();
        });

        this._viewer.addHandler('rotate', function(evt) {
            self.resize();
        });

        this._viewer.addHandler('flip', function() {
          self.resize();
        });

        this._viewer.addHandler('resize', function() {
            self.resize();
        });

        this.resize();
    };

    // ----------
    Overlay.prototype = {
        // ----------
        node: function() {
            return this._node;
        },

        // ----------
        resize: function() {
            if (this._containerWidth !== this._viewer.container.clientWidth) {
                this._containerWidth = this._viewer.container.clientWidth;
                this._svg.setAttribute('width', this._containerWidth);
            }

            if (this._containerHeight !== this._viewer.container.clientHeight) {
                this._containerHeight = this._viewer.container.clientHeight;
                this._svg.setAttribute('height', this._containerHeight);
            }

            var p = this._viewer.viewport.pixelFromPoint(new $.Point(0, 0), true);
            var zoom = this._viewer.viewport.getZoom(true);
            var rotation = this._viewer.viewport.getRotation();
            var flipped = this._viewer.viewport.getFlip();
            // TODO: Expose an accessor for _containerInnerSize in the OSD API so we don't have to use the private variable.
            var containerSizeX = this._viewer.viewport._containerInnerSize.x
            var scaleX = containerSizeX * zoom;
            var scaleY = scaleX;

            if(flipped){
                // Makes the x component of the scale negative to flip the svg
                scaleX = -scaleX;
                // Translates svg back into the correct coordinates when the x scale is made negative.
                p.x = -p.x + containerSizeX;
            }

            this._node.setAttribute('transform',
                'translate(' + p.x + ',' + p.y + ') scale(' + scaleX + ',' + scaleY + ') rotate(' + rotation + ')');
        },
        // ----------
        onClick: function(node, handler) {
            // TODO: Fast click for mobile browsers

            new $.MouseTracker({
                element: node,
                clickHandler: handler
            }).setTracking(true);
        }
    };

})();

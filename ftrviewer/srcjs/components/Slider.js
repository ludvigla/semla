const Slider = (props) => (
  <>
    <input
      onChange={props.handleOpacityChange}
      type="range"
      min={0}
      max={1}
      step={0.01}
      value={props.opacity}
      className="slider"
      name="opacity"
    />
    <label htmlFor="opacity">Opacity</label>
  </>
)

export default Slider

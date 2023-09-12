import { OSG } from './OSG'

const App = (props) => {

  console.log(props)

    return (
      <>
        <OSG
          isNumeric={props.isNumeric}
          width={props.width}
          height={props.height}
          host={props.host}
          port={props.port}
          sampleID={props.sampleID}
          levels={props.levels}
          colors={props.colors}
          scaleByOpacity={props.scaleByOpacity}
          categories={props.categories}
          values={props.values}
          opacities={props.opacities}
          range={props.range}
          useLasso={props.useLasso}
          opacity={props.opacity}
        />
      </>
     )
}

export default App

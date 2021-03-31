####
#
# R Shiny interface to the curatedPCaData-package
#
####

#' Launch R Shiny interface for the curatedPCaData-package
#'
#' Description
#'
#' @rdname shiny
#'
#' @export
shiny <- function(){
	app <- shiny::shinyApp(ui = curatedPCaData:::ui, server = curatedPCaData:::server)
	shiny::runApp(app)
}

#' Shiny UI side for browser
#'
#' Description
#'
#' @rdname shiny
ui <- shiny::navbarPage(
	# Title appended by current package version, sanitized from R
	paste0("curatedPCaData v", utils::packageVersion("curatedPCaData")),
	shiny::tabPanel("Overview",
		shiny::sidebarLayout(
			# Selection menu for exported mae_objects from the package
			shiny::sidebarPanel(
				shiny::uiOutput("maes")
			),
			# Show raw R output for an MAE object
			shiny::mainPanel(
				shiny::verbatimTextOutput("dat_verb")
			)
		)
	),
	shiny::tabPanel("Inspect",
		shiny::sidebarLayout(
			shiny::sidebarPanel(
				"foo"
			),
			shiny::mainPanel(
				"bar"
			)
		)
	)
)

#' Shiny server side
#'
#' Description
#'
#' @rdname shiny
server <- function(input, output, session){
	# Need to load select options on the run
	updateSelectInput(session, "dat",
		choices = as.list(utils::data(package="curatedPCaData")$results[,"Item"])
	)
	# outputs
	output$dat_verb <- shiny::renderPrint({
		if(is.null(input$dat) | input$dat == ""){
			#"Please select a MAE-object exported from the curatedPCaData-package on the left side."
			print(curatedPCaData:::as.named.list(utils::data(package="curatedPCaData")$results[,"Item"]))
		}else{
			eval(parse(text=paste0("curatedPCaData::", input$dat)))
		}
	})
	# MAE set extracted from latest package data listing
	output$maes <- shiny::renderUI({
		shiny::selectInput("dat", "Select MAE-object",
			choices = as.list(utils::data(package="curatedPCaData")$results[,"Item"])
			#choices = list("mae_tcga", "mae_tcga", "mae_sun" = "mae_sun")
		)
	})
}

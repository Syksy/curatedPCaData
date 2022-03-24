####
#
# R Shiny interface to the curatedPCaData-package
#
####

#' Launch R Shiny interface for the curatedPCaData-package
#'
#' Description
#'
#' @rdname shinyPCa
#'
#' @export
shinyPCa <- function(){
	app <- shiny::shinyApp(
		ui = shiny::navbarPage(
			# Title appended by current package version, sanitized from R
			paste0("curatedPCaData v", utils::packageVersion("curatedPCaData")),
			# Show datasets offered in the data set
			shiny::tabPanel("Datasets",
				shiny::sidebarLayout(
					# Selection menu for exported mae_objects from the package
					shiny::sidebarPanel(
						shiny::p("Below are the MultiAssayExperiment (MAE) objects offered in current curatedPCaData package."),
						shiny::p("Select a dataset to further inspect elements contained in the MAE-object."),
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
					# Render a menu for slots available in the currently select MAE-object
					shiny::sidebarPanel(				
						shiny::uiOutput("slots"),
						#shiny::checkboxInput("tSlots", "Transpose", FALSE)
					),
					# Render the selected slot as a DT
					shiny::mainPanel(
						DT::dataTableOutput("slotDT")
					)
				)
			),
			shiny::tabPanel("Clinical",
				shiny::sidebarLayout(
					shiny::sidebarPanel(		
						shiny::p("Clinical data dictionary entries")
						#shiny::uiOutput("slots")
						#shiny::checkboxInput("tClinical", "Transpose", FALSE)
					),
					shiny::mainPanel(
						DT::dataTableOutput("clinicalDT")
					)
				)
			)
		), 
		server = function(input, output, session){
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
				)
			})
			# Data slots available in currently selected data
			output$slots <- shiny::renderUI({
				slots <- eval(parse(text=paste0("names(MultiAssayExperiment::experiments(curatedPCaData::", input$dat,"))")))
				shiny::selectInput("slot", "Select experiment/variables in MAE",
					choices = as.list(slots)
				)		
			})
			# Render a DT (mutable data table with multiple JS features)
			output$slotDT <- DT::renderDataTable({
				#if(input$tSlots){
					eval(parse(text=paste0("as.data.frame(curatedPCaData::", input$dat,"[['", input$slot, "']])")))
				#}else{
				#	eval(parse(text=paste0("as.data.frame(t(curatedPCaData::", input$dat,"[['", input$slot, "']]))")))
				#}
			})
			# Render a DT for clinical variables
			output$clinicalDT <- DT::renderDataTable({
				eval(parse(text=paste0("as.data.frame(as.matrix(MultiAssayExperiment::colData(curatedPCaData::", input$dat,")))")))
			})
		}
	)
	shiny::runApp(app)
}

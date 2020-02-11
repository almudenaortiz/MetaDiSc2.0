
library(shiny)
library(plyr)
library(lmtest)
library(rclipboard)
source("funciones.R")


# START OF SERVER LOGIC
shinyServer(function(input, output, session) {
   
  

  
  rv = reactiveValues()
  
  
  combinations <-  data.frame("Models" = c("A", "B", "C", "D", "E (B')", "F (B')", "G (B')", "H (C')", "I (C')", "J (C')", "K (D')", "L (D')", "M (D')"),
                              "sen" = c("x", "", "", "", "", "", "", "", "", "", "", "", ""),
                              "spe" = c("x", "", "", "", "", "", "", "", "", "", "", "", ""),
                              "sen_0" = c("", "x", "", "x", "x", "x", "x", "", "", "", "x", "x", "x"),
                              "spe_0" = c("", "x", "x", "", "x", "x", "x", "x", "x", "x", "", "", ""),
                              "sen_1" = c("", "x", "", "x", "x", "x", "x", "", "", "", "x", "x", "x"),
                              "spe_1" = c("", "x", "x", "", "x", "x", "x", "x", "x", "x", "", "", ""),
                              "var_sen_0" = c("", "", "", "", "x", "x", "", "x", "x", "", "x", "x", ""),
                              "var_spe_0" = c("", "", "", "", "x", "", "x", "x", "", "x", "x", "", "x"),
                              "var_sen_1" = c("", "", "", "", "x", "x", "", "x", "x", "", "x", "x", ""),
                              "var_spe_1" = c("", "", "", "", "x", "", "x", "x", "", "x", "x", "", "x")
  )
  
  
  output$combinations <- renderTable(
    combinations, align = 'c', striped = TRUE, bordered = T, spacing = "s")
  
    
  
  # CSV FILE AS INPUT
  observeEvent(input$file1, ({
    rv$datos = read.csv(input$file1$datapath)
   
    
    colnames(rv$datos) = tolower(colnames(rv$datos))
    
    if (!c("id", "tp", "tn", "fn", "fp") %in% colnames(rv$datos)){
      showModal(modalDialog("Required columns are missing"))
      return()
    }
    
    if (!grep("id", colnames(rv$datos)))
    {
      showNotification("Error: The data does not contain a 'ID' column", type = "error")
      return()
    } else if (!grep("tp", colnames(rv$datos)))
    {
      showNotification("Error: The data does not contain a 'TP' column", type = "error")
      return()
    } else if (!grep("tn", colnames(rv$datos)))
    {
      showNotification("Error: The data does not contain a 'TN' column", type = "error")
      return()
    } else if (!grep("fp", colnames(rv$datos)))
    {
      showNotification("Error: The data does not contain a 'FP' column", type = "error")
      return()
    } else if (!grep("fn", colnames(rv$datos)))
    {
      showNotification("Error: The data does not contain a 'FN' column", type = "error")
      return()
    } 
    
    output$DataImport = DT::renderDataTable({
      datatable(rv$datos, rownames = F
                # ,
                # extensions = c("Buttons", "Select"),
                # options = list(
                #   select = TRUE,
                #   dom = "Bfrtip",
                #   buttons = list(
                #     list(
                #       extend = "copy",
                #       text = "Copy",
                #       exportOptions = list(modifier = list(selected=TRUE))
                #     )
                #   )
                # )
                )
    }, server = TRUE)
    

    
  }))
  
  
  
  observeEvent(input$Apply, ({
    if (is.null(rv$datos))
    {
      showNotification("Error: Not file loaded", type = "error")
      return()
    } 
    
    
    withProgress(message = "Performing Analysis", {
      
      datos <- rv$datos[input$DataImport_rows_all, ]
      

      
      #selection of covariates
      select_subgroup <- reactiveValues() 
      select_subgroup <- colnames(datos %>% select(-c("id","tp","fp","tn","fn")))
      updateSelectInput(session, "subgroups_list", choices = select_subgroup)
      
      updateSelectInput(session, "covariates_list", choices = select_subgroup)
      
      X = rv$datos[input$DataImport_rows_all, ] %>% mutate(
          n1 = tp + fn,
          n0 = fp + tn,
          true1 = tp,
          true0 = tn,
          study = 1:nrow(rv$datos)
      )
      
      colnames(X) <- tolower(colnames(X))
      X$id <-  make.names(X$id, unique = T, allow_ = F)
      
      
      
      #selection of studies
      selection_studies <-  reactiveValues() 
      selection_studies <-  X$id
      updateCheckboxGroupInput(session, "studies_list", choices = selection_studies, selected = c(selection_studies))
      
      
      ## STUDY-LEVEL OUTCOME
      ##
      
      X2 <- rv$datos[input$DataImport_rows_all, ] %>% mutate(
        n1 = tp + fn,
        n0 = fp + tn,
        pos = tp + fp,
        neg = fn + tn,
        sens = round(tp/n1,3),
        spec = round(tn/n0,3)
      )
      
      output$slevelTable <- DT::renderDataTable(datatable(X2 %>% select(id, tp, fp, fn, tn, n1, n0, pos, neg, sens, spec)))
      
      output$downloadData1 <-  downloadHandler(
        filename = function(){
          paste('study-level', Sys.Date(), '.csv', sep = '')
        },
        content = function(con){
          write.csv(X2, con)
        }
      )
      

      
      
      ###
      ### CHANGE DATA
      ###
      
      
      Y <- formato_datos(datos)
      
      
      ###
      ### FOREST PLOTS of sens and spec
      ###
      
      
      Y_forest <-  format_forest(X)

      output$CIsen = renderPlotly({
        forest_sens(Y_forest)
      })
      
      # output$downloadCIsen <- downloadHandler(
      #   filename = "forest_plot_sens.png",
      #   content = function(file){
      #     file.copy(forest_sens(Y_forest), file)
      #   }
      # )
      
      output$CIspe = renderPlotly({
        forest_spec(Y_forest)
      })
      
      
      
      
      
      ###
      ### BIVARIATE BINOMIAL META-ANALYSIS WITH GLMER
      ### 
      
      
      ## GLM table estimates
      
      ma_biv = modelo(X)$summary
      logit <- modelo(X)$logit
      tabla <- modelo(X)$tabla
      corr <- modelo(X)$corr
      

      output$glm_table <- renderTable(
        tabla,  striped = TRUE, digits = 3
      )
      
      output$downloadDataGlm <-  downloadHandler(
        filename = function(){
          paste('glm_estimates', Sys.Date(), '.csv', sep = '')
        },
        content = function(con){
          write.csv(tabla, con)
        }
      )
      
      
      ## RevMan parameters table
      
      corr_df <- data.frame(
        Coefficient = "Corr(logits)",
        Estimate = corr
      )
      
      tablarevman <- plot_sroc(Y, ma_biv)$tabla
      tablarevman2 <- rbind(tablarevman, corr_df)
      
      output$tablaRevMan <- renderTable(
        tablarevman2,  striped = TRUE, digits = 3
      )
      
      output$downloadRevMan <-  downloadHandler(
        filename = function(){
          paste('revman_parameters', Sys.Date(), '.csv', sep = '')
        },
        content = function(con){
          write.csv(tablarevman2, con)
        }
      )
      

      
      tabla_ellipse <- plot_sroc(Y, ma_biv)$tabla_ellipse
      
      
      
      ###
      ### SROC plot options
      ###
      
      
      observeEvent(input$sroc_options, {
        sroc_plot <- plot_ly() %>% layout(
          title = "SROC plane",
          xaxis = list(title = "1-Specificity", range = c(-0.05, 1.05)),
          yaxis = list(title = "Sensitivity", range = c(-0.05, 1.05)),
          showlegend = FALSE,
          width = 450, 
          height = 450
        )
        if (identical(input$sroc_options, "studies")){
                output$SROC <- renderPlotly({
                  sroc_plot %>% add_trace(
                    data = Y %>% mutate(x = 1 - (tn / (tn + fp)), y = tp / (tp + fn)),
                    x =  ~ x,
                    y =  ~ y,
                    type = "scatter",
                    hoverinfo = "text",
                    text = ~ paste('ID:', id, '\n', '1-Specificity:', x, '\n', 'Sensitivity:' , y),
                    marker = list(color = pal_npg("nrc")(1))
                  )
                })
        } else if (identical(input$sroc_options, "ellipse")) {
              output$SROC <- renderPlotly({
                sroc_plot %>% add_trace(
                  data = tabla_ellipse,
                  x = ~ specificity,
                  y = ~ sensitivity,
                  text = NULL,
                  mode = "lines",
                  marker = list(size = 0.5, color = toRGB("black"))
                )
              })
        } else {
                output$SROC <- renderPlotly({
                  sroc_plot = plot_sroc(Y, ma_biv)$sroc
                  sroc_plot
                })
        }
        
      })
      

      
      
      
      ###
      ### HETEROGENEITY ESTIMATES TABLE BY BAYESIAN MODEL
      ###
      


      area <- round(plot_sroc(Y, ma_biv)$area,3)
      VarLogitSen = round(ma_biv$varcor$study[1, 1],3)  ##Var(logit(sen))
      VarLogitSpe = round(ma_biv$varcor$study[2, 2],3)  ##Var(logit(spe))
      
      tabla_het2 = data.frame(
        Coefficient = c("Area of ellipse", "Var(logit(sen))", "Var(logit(spec))"),
        Estimate = c(area, VarLogitSen, VarLogitSpe)
      )
      
      
      output$tabla_het <-  renderTable(
        tabla_het2, align = 'c', striped = TRUE, digits = 3
      )
      
      
      output$downloadHet_est <-  downloadHandler(
        filename = function(){
          paste('heterogeneity_estimates', Sys.Date(), '.csv', sep = '')
        },
        content = function(con){
          write.csv(tabla_het2, con)
        }
      )
      
      
      observeEvent(input$heter, {
        withProgress(message = "Calculating heterogeneity estimates", detail = "This may take a while...", {
          heterogen <- heterogeneidad(X)
          
          output$Table2 = renderTable(
            heterogen, align = 'c', striped = TRUE, digits = 3)
          
          output$downloadDataHet <-  downloadHandler(
            filename = function(){
              paste('heterogeneity_estimates', Sys.Date(), '.csv', sep = '')
            },
            content = function(con){
              write.csv(heterogen, con)
            }
          )
        })
        })
      
      
      
      
      
      ###
      ### ANALYSIS OF SUBGROUPS
      ###
  
      
      
      observeEvent(input$subgroup_button, {
        subgroup1  <- paste("Subgroup selected: ",as.character(input$subgroups_list))
        subgroup <-  input$subgroups_list
        output$subgroup <- renderText(subgroup1)
        
        X3 <- split(X, as.factor(X[,paste(subgroup)]))
        modelo_covariables <- modelo_sg(X3)
        names(modelo_covariables) = names(X3)
        
        
        ## Estimates Table

        output$glm_subgroups <- renderUI({
          tfc = function(m, name){
            renderDataTable(m, options = list(paging = FALSE, searching = FALSE), caption = paste("Subgroup: ",name,"."))
          }
          list_of_tables <-  lapply(1:length(modelo_covariables), function(i){
            name_subgroup = names(modelo_covariables)[i]
            tfc(modelo_covariables[[i]]$resultados_final$df, name_subgroup)
          })
        })
        

        
        
        ## SROC plots
        
        df <-  ldply(X3, data.frame)
        df$Specif <- 1-(df$tn/(df$tn+df$fp))
        df$Sensit <- df$tp/(df$tp+df$fn)
        
        plot_tmp <-  ggplot() + geom_point(data = df, aes(x= Specif, y = Sensit, text= paste("ID:", id),  colour = factor(df[,(subgroup)]))) 
        
        mypal <-  pal_npg("nrc")(length(modelo_covariables))
        
        
        for (i in 1:length(modelo_covariables)) {
          
          plot_tmp <- plot_tmp   + 
            geom_polygon(data = modelo_covariables[[i]]$tabla2$tabla2, aes(specificity, sensitivity),fill = NA, colour = mypal[i]) + scale_color_npg() 
        }
        
        plot_tmp <-   plot_tmp + xlim(0,1) + ylim(0,1) +  
          labs(x = "1-Specificity", y = "Sensitivity", title= "SROC Plane", colour ="Covariates") 
        
        
        
        output$sroc_subgroups <- renderPlotly({
          ggplotly(plot_tmp, tooltip = c("x","y","text"))
        }) 
         
        
      })
      
      
      
      
      
      
      
      ###
      ### META-REGRESSION
      ###
      
      ## DIfferent models table 

      
      ## Meta-reg button 
      observeEvent(input$metareg_button, {
        
        
        selected = input$covariates_list
        
        
        X_model <- X
        
        X_model$mr <- X[,paste(selected)]
        

        if (nlevels(X_model$mr) == 2) {
        
          nsub_texto = paste("Subgroup 0 = ", levels(X[,paste(selected)])[1], ", Subgroup 1 = ", levels(X[,paste(selected)])[2])
          
          output$nsub = renderPrint({
            cat(nsub_texto)
          })
          
          
          X_model$mr <- as.numeric(X_model$mr)
          X_model$mr <- (X_model$mr)-1
          X_model <-  reshape(X_model, direction = "long", varying = list(c("n1", "n0"), c("true1", "true0")), timevar = "sens", times = c(1,0),
                              v.names = c("n", "true"))
          
          X_model <- X_model[order(X_model$id),]
          X_model$spec <-  1-X_model$sens
          X_model$A <- 1-(X_model$mr)
          X_model$B <- 1-X_model$A
          
          X_model$sen_0 <- (X_model$A)*(X_model$sens)
          X_model$sen_1 <- (X_model$B)*(X_model$sens)
          X_model$spe_0 <- (X_model$A)*(X_model$spec)
          X_model$spe_1 <- (X_model$B)*(X_model$spec)
          
          ## CASCADA MODELOS
          
          convg_A <- function(X_model){ # model A
            out <- tryCatch(
              {
                model_A = glmer(formula = cbind(true, n - true)~ 0 + sens + spec + (0+sens + spec|study), data = X_model, family = binomial)
              }, 
              error=function(cond){
                message("Here is the original error message:")
                message(cond)
                return(NA)
              },
              warning=function(cond){
                converge = "Model A failed to converge."
                message("Here is the original warning message:")
                message(cond)
                return(NULL)
              }
              
            )
            return(out)
          } # end model A
          
          model_A = convg_A(X_model)
          
          convg_B <- function(X_model){ # model B 
            out <- tryCatch(
              {
                model_B = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spe_0 + spe_1 + (0+sens + spec|study), data = X_model, family = binomial)
              }, 
              error=function(cond){
                message("Here is the original error message:")
                message(cond)
                return(NA)
              },
              warning=function(cond){
                converge = "Model B failed to converge."
                message("Here is the original warning message:")
                message(cond)
                return(NULL)
              }
            )
            return(out)
          } # end model B 
          
          model_B = convg_B(X_model)
          test_BA = lrtest(model_A, model_B)
          pval_BA <- test_BA$`Pr(>Chisq)`[2]
          
          final = ""
          
          if (pval_BA <= 0.05) {  # si signif AB
            
            # model_C
            convg_C <- function(X_model){
              out <- tryCatch(
                {
                  model_C = glmer(formula = cbind(true, n - true)~ 0 + sens + spe_0 + spe_1 + (0+sens + spec|study), data = X_model, family = binomial)
                }, 
                error=function(cond){
                  message("Here is the original error message:")
                  message(cond)
                  return(NA)
                },
                warning=function(cond){
                  converge="Model C failed to converge."
                  message("Here is the original warning message:")
                  message(cond)
                  return(NULL)
                }
                
              )
              return(out)
            } # end model C
            model_C = convg_C(X_model)
            test_BC = lrtest(model_B, model_C)
            pval_BC <- test_BC$`Pr(>Chisq)`[2]
            
            #model_D
            
            convg_D <- function(X_model){
              out <- tryCatch(
                {
                  model_D = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spec + (0+sens + spec|study), data = X_model, family = binomial)
                }, 
                error=function(cond){
                  message("Here is the original error message:")
                  message(cond)
                  return(NA)
                },
                warning=function(cond){
                  converge="Model D failed to converge."
                  message("Here is the original warning message:")
                  message(cond)
                  return(NULL)
                }
                
              )
              return(out)
            } # end model D
            model_D = convg_D(X_model)
            test_BD = lrtest(model_B, model_D)
            pval_BD <- test_BD$`Pr(>Chisq)`[2]
            
            if (pval_BC <= 0.05 & pval_BD > 0.05) { #si signif BC, no signif BD - model_D
              
              # model_K
              convg_K <- function(X_model){
                out <- tryCatch(
                  {
                    model_K = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spec +(0 +sen_1 + spe_1 |study) + (0 +sen_0 + spe_0 |study), data = X_model, family = binomial)
                  }, 
                  error=function(cond){
                    message("Here is the original error message:")
                    message(cond)
                    return(NA)
                  },
                  warning=function(cond){
                    converge="Model K failed to converge."
                    message("Here is the original warning message:")
                    message(cond)
                    return(NULL)
                  }
                  
                )
                return(out)
              } # end model K 
              model_K = convg_K(X_model)
              test_DK = lrtest(model_D, model_K)
              pval_DK <- test_DK$`Pr(>Chisq)`[2]
              
              if (pval_DK <= 0.05) { #si signif DK
                
                # model_L
                convg_L <- function(X_model){
                  out <- tryCatch(
                    {
                      model_L = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spec +(0 +sen_1 + spec |study) + (0 +sen_0 + spec |study), data = X_model, family = binomial)
                    }, 
                    error=function(cond){
                      message("Here is the original error message:")
                      message(cond)
                      return(NA)
                    },
                    warning=function(cond){
                      converge="Model L failed to converge."
                      message("Here is the original warning message:")
                      message(cond)
                      return(NULL)
                    }
                    
                  )
                  return(out)
                } # end model L 
                model_L = convg_L(X_model)
                test_KL = lrtest(model_K, model_L)
                pval_KL <- test_KL$`Pr(>Chisq)`[2]
                
                # model_M
                
                convg_M <- function(X_model){
                  out <- tryCatch(
                    {
                      model_M = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spec +(0 +spe_1 + sens |study) + (0 +sens + spe_0 |study), data = X_model, family = binomial)
                    }, 
                    error=function(cond){
                      message("Here is the original error message:")
                      message(cond)
                      return(NA)
                    },
                    warning=function(cond){
                      converge="Model M failed to converge."
                      message("Here is the original warning message:")
                      message(cond)
                      return(NULL)
                    }
                    
                  )
                  return(out)
                } # end model M 
                model_M = convg_M(X_model)
                test_KM = lrtest(model_K, model_M)
                pval_KM <- test_KM$`Pr(>Chisq)`[2]
                
                
                if (pval_KL <= 0.05 & pval_KM > 0.05) { #si signif KL, no signif KM
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and C but not between B and D.\nThere is statistical significance between models D and K.\nThere is statistical significance between models K and L but not between models K and M.\nThe final model is M."
                  model_final = model_M
                } else if (pval_KM <= 0.05 & pval_KL > 0.05) { #si signif KM, no signif KL
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and C but not between B and D.\nThere is statistical significance between models D and K.\nThere is statistical significance between models K and M but not between models K and L.\nThe final model is L."
                  model_final = model_L
                } else if (pval_KL <= 0.05 & pval_KM <= 0.05) { #si signif KM, si signif KL
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and C but not between B and D.\nThere is statistical significance between models D and K.\nThere is statistical significance between models K and L, and also between models K and M.\nThe final model is K."
                  model_final = model_K
                }
                
              } else { # no signif DK
                
                if (is.null(model_K)){
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and D (pVal < 0.05).\nModel K failed to converge.\nThe final model is D."
                }else {
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and D (pVal < 0.05).\nThere is NOT statistical significance between models D and K.\nThe final model is D."
                }
                
                model_final = model_D
              }
              
            } else if (pval_BD <= 0.05 & pval_BC > 0.05) { #si signif BD, no signif BC - model_C
              
              # model_H
              convg_H <- function(X_model){
                out <- tryCatch(
                  {
                    model_H = glmer(formula = cbind(true, n - true)~ 0 + sens + spe_0 + spe_1 +(0 +sen_1 + spe_1 |study) + (0 +sen_0 + spe_0 |study), data = X_model, family = binomial)
                  }, 
                  error=function(cond){
                    message("Here is the original error message:")
                    message(cond)
                    return(NA)
                  },
                  warning=function(cond){
                    converge="Model H failed to converge."
                    message("Here is the original warning message:")
                    message(cond)
                    return(NULL)
                  }
                  
                )
                return(out)
              } # end model H 
              model_H = convg_H(X_model)
              test_CH = lrtest(model_C, model_H)
              pval_CH <- test_CH$`Pr(>Chisq)`[2]
              
              
              if (pval_CH <= 0.05){ #si signif CH
                
                # model_I
                convg_I <- function(X_model){
                  out <- tryCatch(
                    {
                      model_I = glmer(formula = cbind(true, n - true)~ 0 + sens + spe_0 + spe_1 +(0 +sen_1 + spec |study) + (0 +sen_0 + spec |study), data = X_model, family = binomial)
                    }, 
                    error=function(cond){
                      message("Here is the original error message:")
                      message(cond)
                      return(NA)
                    },
                    warning=function(cond){
                      converge="Model I failed to converge."
                      message("Here is the original warning message:")
                      message(cond)
                      return(NULL)
                    }
                    
                  )
                  return(out)
                }  # end model I 
                model_I = convg_I(X_model)
                test_HI = lrtest(model_H, model_I)
                pval_HI <- test_HI$`Pr(>Chisq)`[2]
                
                # model_J
                convg_J <- function(X_model){
                  out <- tryCatch(
                    {
                      model_J = glmer(formula = cbind(true, n - true)~ 0 + sens + spe_0 + spe_1 +(0 +spe_1 + sens |study) + (0 +sens + spe_0 |study), data = X_model, family = binomial)
                    }, 
                    error=function(cond){
                      message("Here is the original error message:")
                      message(cond)
                      return(NA)
                    },
                    warning=function(cond){
                      converge="Model J failed to converge."
                      message("Here is the original warning message:")
                      message(cond)
                      return(NULL)
                    }
                    
                  )
                  return(out)
                } # end model J 
                model_J = convg_J(X_model)
                test_HJ = lrtest(model_H, model_J)
                pval_HJ <- test_HJ$`Pr(>Chisq)`[2]
                
                
                if (pval_HI <= 0.05 & pval_HJ > 0.05) { #si signif HI, no signif HJ
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and D but not between B and C.\nThere is statistical significance between models C and H.\nThere is statistical significance between models H and I, but not between models H and J.\nThe final model is J."
                  model_final = model_J
                } else if (pval_HJ <= 0.05 & pval_HI > 0.05) { #si signif HJ, no signif HI
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and D but not between B and C.\nThere is statistical significance between models C and H.\nThere is statistical significance between models H and J, but not between models H and I.\nThe final model is I."
                  model_final = model_I
                } else if(pval_HI <= 0.05 & pval_HJ <= 0.05) { #si signif HI, si signif HJ
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and D but not between B and C.\nThere is statistical significance between models C and H.\nThere is statistical significance between models H and I, and also between models H and J.\nThe final model is H."
                  model_final = model_H
                }
                
              } else { #no signif CH
                if (is.null(model_H)) {
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and D but not between B and C.\nModel H failed to converge.\nThe final model is C."
                } else {
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is statistical significance between models B and D but not between B and C.\nThere is NOT statistical significance between models C and H.\nThe final model is C."
                }
                
                model_final = model_C
              }
              
            } else if (pval_BD <= 0.05 & pval_BC <= 0.05){ #si signif BD y si signif BC - model_B
              
              # model_E
              convg_E <- function(X_model){
                out <- tryCatch(
                  {
                    model_E = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spe_0 + spe_1 +(0 +sen_1 + spe_1 |study) + (0 +sen_0 + spe_0 |study), data = X_model, family = binomial)
                  }, 
                  error=function(cond){
                    message("Here is the original error message:")
                    message(cond)
                    return(NA)
                  },
                  warning=function(cond){
                    converge = "Model E failed to converge."
                    message("Here is the original warning message:")
                    message(cond)
                    return(NULL)
                  }
                  
                )
                return(out)
              } # end model E
              model_E = convg_E(X_model)
              test_EB = lrtest(model_B, model_E)
              pval_EB <- test_EB$`Pr(>Chisq)`[2]
              
              if (pval_EB <= 0.05) { #si signif EB
                
                # model_F
                convg_F <- function(X_model){
                  out <- tryCatch(
                    {
                      model_F = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spe_0 + spe_1 +(0 +sen_1 + spec |study) + (0 +sen_0 + spec |study), data = X_model, family = binomial)
                    }, 
                    error=function(cond){
                      message("Here is the original error message:")
                      message(cond)
                      return(NA)
                    },
                    warning=function(cond){
                      converge="Model F failed to converge."
                      message("Here is the original warning message:")
                      message(cond)
                      return(NULL)
                    }
                    
                  )
                  return(out)
                } # end model F 
                model_F = convg_F(X_model)
                test_EF = lrtest(model_E, model_F)
                pval_EF <- test_EF$`Pr(>Chisq)`[2]
                
                # model_G
                convg_G <- function(X_model){
                  out <- tryCatch(
                    {
                      model_G = glmer(formula = cbind(true, n - true)~ 0 + sen_0 + sen_1 + spe_0 + spe_1 +(0 +spe_1 + sens |study) + (0 +spe_0 + sens |study), data = X_model, family = binomial)
                    }, 
                    error=function(cond){
                      message("Here is the original error message:")
                      message(cond)
                      return(NA)
                    },
                    warning=function(cond){
                      converge="Model G failed to converge."
                      message("Here is the original warning message:")
                      message(cond)
                      return(NULL)
                    }
                    
                  )
                  return(out)
                } # end model G 
                model_G = convg_G(X_model)
                test_EG = lrtest(model_E, model_G)
                pval_EG <- test_EG$`Pr(>Chisq)`[2]
                
                
                if(pval_EF <= 0.05 & pval_EG > 0.05){ #si signif EF, no signific EG
                  final = "There is statistical significance between models A and B (pVal < 0.05).\n There is statistical significance between models B and C, and also B and D.\nThere is statistical difference between models E and B.\nThere is statistical difference between models E and F, but not between E and G.\n The final model is G. "
                  model_final = model_G
                } else if (pval_EG <= 0.05 & pval_EF > 0.05){ #si signif EG, no signif EF
                  final = "There is statistical significance between models A and B (pVal < 0.05).\n There is statistical significance between models B and C, and also B and D.\nThere is statistical difference between models E and B.\nThere is statistical difference between models E and G, but not between E and F.\n The final model is F. "
                  model_final = model_F
                } else if (pval_EF <= 0.05 & pval_EG <= 0.05){ #si signif EF, si signif EG
                  final = "There is statistical significance between models A and B (pVal < 0.05).\n There is statistical significance between models B and C, and also B and D.\nThere is statistical difference between models E and B.\nThere is statistical difference between models E and F, and also E and G.\n The final model is E. "
                  model_final = model_E
                }
                
              }  else { # no signif EB
                
                if (is.null(model_E)){
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nModel E failed to converge.\nThe final model is B."
                } else {
                  final = "There is statistical significance between models A and B (pVal < 0.05).\nThere is NOT statistical significance between models B and E (pVal > 0.05).\nThe final model is B."
                }
                
                #final = "model_B"
                model_final = model_B
              }
              
            } 
            
          } else { # no signif AB
            
            if (is.null(model_B)){
              final = "There is statistical significance between models A and B (pVal < 0.05).\nModel B failed to converge.\nThe final model is A."
            } else {
              final = "There is NOT statistical significance between models A and B (pVal > 0.05).\nThe final model is A."
            }
            
            model_final = model_A
          }
          
          ## fin cascada
          
          output$models <- renderPrint({
            cat(final)
          })
          prueba = mr_table(model_A, model_final)$tabla
          
          output$model_final <- renderTable({
            prueba
          })
          
        } else {
          showModal(modalDialog("Select a covariate with two factor levels."))
          return()
        }
        

        
        
        
      }) # fin button meta-reg
      
      
      
        
      ### 
      ### SENSIBILITY
      ### 
      
      observeEvent(input$studies_list, {
        
        # Selected studies
        
        X2_studies <- X[which(X$id %in% input$studies_list),] #table with information of selected studies only
        
        
        
        observeEvent(input$refresh, {
          ## GLM table estimates
          ma_biv_sens = modelo(X2_studies)$summary
          logit_sens <- modelo(X2_studies)$logit
          tabla_sens <- modelo(X2_studies)$tabla
          corr_sens <- modelo(X2_studies)$corr
          
          output$glm_table_sens <- renderTable(
            tabla_sens,  striped = TRUE, digits = 3
          )
          
          output$downloadDataGlm_sens <-  downloadHandler(
            filename = function(){
              paste('glm_estimates', Sys.Date(), '.csv', sep = '')
            },
            content = function(con){
              write.csv(tabla_sens, con)
            }
          )
          
          ## SROC Plane
          
          Y <- formato_datos(X2_studies)
          
          sroc_sens <- plot_sroc(Y, ma_biv_sens)$sroc
          
          output$sroc_sens <- renderPlotly({
            sroc_sens
          })
          
        })
        
      })
      
      
      


      updateTabsetPanel(session, "tabSet", selected = "Plots")
    } )
  }))
  
})

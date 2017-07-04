  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  module m_dag_interpol
    use m_dag_interpol_0
    use m_dag_interpol_1
    use m_dag_interpol_2
    use m_dag_interpol_3
    use m_dag_interpol_dim
    use m_dag_interpol_dim_2
    implicit none
    
    interface interpol_1o
      module procedure interpol_1o, interpol_1o_r
    end interface
    
    interface interpol_2o
      module procedure interpol_2o, interpol_2o_r
    end interface
    
    interface interpol_linear_1D
      module procedure interpol_linear_1D, interpol_linear_1D_r, &
                       interpol_linear_1D_2, interpol_linear_1D_r_2
    end interface
    
    interface interpol_quadratic_1D
      module procedure interpol_quadratic_1D, interpol_quadratic_1D_r, &
                       interpol_quadratic_1D_2, interpol_quadratic_1D_r_2
    end interface
    
    interface interpol_linear_2D
      module procedure interpol_linear_2D, interpol_linear_2D_r, &
                       interpol_linear_2D_3, interpol_linear_2D_r_3
    end interface
    
    interface interpol_linear_3D
      module procedure interpol_linear_3D, interpol_linear_3D_r
    end interface
    
  !***********************************************************************/
  end module m_dag_interpol
  
  
  
  
  
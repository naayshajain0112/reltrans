! Import FTorch library for interfacing with PyTorch
use ftorch, only : torch_model, torch_tensor, torch_kCPU, torch_delete, &
                   torch_tensor_from_array, torch_model_load, torch_model_forward

implicit none

subroutine init_emulator()

    ! Generate an object to hold the Torch model
    type(torch_model) :: emulator
    
    ! Initialise the Torch model to be used
    torch_model_load(emulator, pathname_emulator, torch_kCPU)

end subroutine init_emulator

subroutine run_model_inference(emulator, photar, nex)

    ! Set up array of n_inputs input tensors and array of n_outputs output tensors
    ! Note: In this example there is only one input tensor (n_inputs = 1) and one
    !       output tensor (n_outputs = 1)
    integer, parameter :: n_inputs = 10
    integer, parameter :: n_outputs = 1
    type(torch_tensor), dimension(n_inputs)  :: model_pars
    type(torch_tensor), dimension(n_inputs)  :: model_egrid
    type(torch_tensor), dimension(n_outputs) :: model_output_arr

    ! Set up the model inputs and output as Fortran arrays
    real, dimension(10), target    :: pars
    real, dimension(nex), target   :: egrid
    real, dimension(nex), target     :: output

    ! Initialise the inputs as Fortran array of ones
    input = pars
    do i = 1, nex
        egrid(i) = ear(i)
    end do

    ! Wrap Fortran data as no-copy Torch Tensors
    call torch_tensor_from_array(model_pars(1), input, torch_kCPU)
    call torch_tensor_from_array(model_egrid(1), egrid, torch_kCPU)
    call torch_tensor_from_array(model_output_arr(1), output, torch_kCPU)

    ! Run model forward method and infer reflection spectrum
    call torch_model_forward(emulator, model_pars, model_egrid, model_output_arr)

    ! Clean up
    call torch_delete(model_pars)
    call torch_delete(model_egrid)
    call torch_delete(model_output_arr)

    photar = output

end subroutine run_model_inference
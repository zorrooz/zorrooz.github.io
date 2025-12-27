<!-- AppFooter.vue -->
<template>
  <footer class="py-4 site-footer" style="background-color: var(--app-footer-bg);">
    <div class=" container">
      <div class="d-flex justify-content-between align-items-center">
        <div>Copyright © 2025 gblog</div>
        <div class="d-flex flex-column align-items-center">
          <div> Designed by <span class="fw-bold">zorrooz</span></div>
          <div class="d-flex justify-content-around">
            <span>
              <a href="https://github.com/zorrooz" class="me-3" style="color: var(--app-text-muted);">
                <i class="fab fa-github fa-lg"></i>
              </a>
            </span>
            <span>
              <a href="mailto:zorrooz@163.com" style="color: var(--app-text-muted);">
                <i class="fas fa-envelope fa-lg"></i>
              </a>
            </span>
          </div>
        </div>
      </div>
    </div>
  </footer>
</template>

<script>
/*
  AppFooter
  - 页面底部信息
*/
export default {
  name: 'AppFooter'
}
</script>

<style scoped></style>

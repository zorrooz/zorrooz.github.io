<!-- ProfileCard.vue -->
<template>
  <div class="card shadow-sm border-0 bg-white" style="max-width: 280px;">
    <div class="text-center mt-3">
      <!-- 头像占位 -->
      <div class="bg-light rounded-circle d-inline-flex align-items-center justify-content-center"
        style="width: 80px; height: 80px;">
        <span class="text-muted">头像</span>
      </div>
    </div>

    <div class="card-body p-3 text-center">
      <h6 class="card-title mb-1 fw-bold">zorrooz</h6>
      <p class="card-text text-muted small mb-3">
        开发者
      </p>

      <div class="row g-0 text-center">
        <div class="col border-end">
          <div class="fw-bold">128</div>
          <div class="text-muted small">文章</div>
        </div>
        <div class="col border-end">
          <div class="fw-bold">8</div>
          <div class="text-muted small">分类</div>
        </div>
        <div class="col border-end">
          <div class="fw-bold">42</div>
          <div class="text-muted small">标签</div>
        </div>
        <div class="col">
          <div class="fw-bold">15.6万</div>
          <div class="text-muted small">字数</div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'UserInfoCard'
}
</script>

<style scoped>
/* 无自定义CSS */
</style>

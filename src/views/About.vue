<!-- About.vue -->
<template>
  <div class="container py-4 px-3">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8 typography-body">

        <!-- 页面标题 -->
        <div class="text-center mb-4">
          <h1 class="article-title text-primary mb-3">{{ pageTitle }}</h1>
        </div>

        <!-- 内容区块 -->
        <div class="d-flex flex-column" style="gap: 2rem;">

          <!-- 自我介绍 -->
          <section>
            <h2 class="h4 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              自我介绍
            </h2>
            <p class="text-muted mb-0">
              {{ data.introduction }}
            </p>
          </section>

          <!-- 通用分组（统一 section 列表） -->
          <section v-for="(entry, idx) in data.section" :key="idx">
            <h2 class="h4 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              {{ entry.title }}
            </h2>
            <div class="ms-3">
              <div v-for="(it, j) in entry.items" :key="j" class="mb-4">
                <p class="text-muted fw-bold mb-1">{{ it.item }}</p>
                <p v-if="it.desc" class="text-muted mb-0">{{ it.desc }}</p>
              </div>
            </div>
          </section>

          <!-- 联系我 -->
          <section>
            <h2 class="h4 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              联系我
            </h2>
            <div class="text-center">
              <div class="d-flex flex-wrap justify-content-center align-items-center gap-4 text-muted">
                <a v-for="contact in data.contacts" :key="contact.label" :href="contact.link" target="_blank"
                  rel="noopener noreferrer"
                  class="d-flex align-items-center text-decoration-none text-secondary hover-primary"
                  @mouseenter="e => e.target.style.transform = 'translateY(-1px)'"
                  @mouseleave="e => e.target.style.transform = 'translateY(0)'"
                  :style="{ transition: 'transform 0.2s ease', color: '#6c757d' }">
                  <i :class="contact.icon" class="me-1"></i>
                  <span>{{ contact.value }}</span>
                </a>
              </div>
            </div>
          </section>

        </div>

        <!-- 结尾标语 -->
        <div class="text-center mt-5 pt-4" style="border-top: 1px solid #dee2e6;">
          <p class="text-muted mb-0">
            <i class="fas fa-star me-1"></i>
            {{ footerText }}
          </p>
        </div>

      </div>
    </div>
  </div>
</template>

<script>
import aboutData from '@/content/about.json'

export default {
  name: 'AboutView',
  data() {
    return {
      data: aboutData,
      pageTitle: '关于我',
      footerText: '感谢您的关注！'
    }
  }
}
</script>

<style scoped>
.typography-body { font-size: 1.125rem; line-height: 1.8; color: var(--bs-gray-800); }
.typography-body p { margin-bottom: 0.75rem; }
</style>

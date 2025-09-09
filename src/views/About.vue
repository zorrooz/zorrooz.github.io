<!-- About.vue -->
<template>
  <div class="container py-4">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8">

        <!-- 页面标题 -->
        <div class="text-center mb-4">
          <h1 class="display-6 fw-bold text-primary">{{ data.pageTitle }}</h1>
        </div>

        <!-- 内容区块 -->
        <div class="d-flex flex-column" style="gap: 2rem;">

          <!-- 自我介绍 -->
          <section>
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              自我介绍
            </h2>
            <p class="text-muted mb-0" style="font-size: 1rem; line-height: 1.6; color: #6c757d;">
              {{ data.introduction }}
            </p>
          </section>

          <!-- 教育经历 -->
          <section>
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              教育经历
            </h2>
            <div class="ms-3">
              <div v-for="(edu, idx) in data.education" :key="idx" class="mb-4">
                <h3 class="h6 fw-semibold text-dark mb-1">{{ edu.degree }}</h3>
                <p class="text-muted mb-1 small">
                  <strong>{{ edu.institution }}</strong> · {{ edu.duration }}
                </p>
                <p class="text-muted mb-0 small" style="color: #6c757d;">{{ edu.desc }}</p>
              </div>
            </div>
          </section>

          <!-- 专业技能 -->
          <section>
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              专业技能
            </h2>
            <div class="ms-3">
              <div v-for="(skill, idx) in data.skills" :key="idx" class="mb-4">
                <h3 class="h6 fw-semibold text-dark mb-1">{{ skill.name }}</h3>
                <p class="text-muted mb-1 small">{{ skill.level }}</p>
                <p class="text-muted mb-0 small" style="color: #6c757d;">{{ skill.tools }}</p>
              </div>
            </div>
          </section>

          <!-- 项目经历 -->
          <section>
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              项目经历
            </h2>
            <div class="ms-3">
              <div v-for="(project, idx) in data.projects" :key="idx" class="mb-4">
                <h3 class="h6 fw-semibold text-dark mb-1">{{ project.title }}</h3>
                <p class="text-muted mb-1 small">
                  <strong>{{ project.role }}</strong> · {{ project.duration }}
                </p>
                <p class="text-muted mb-0 small" style="color: #6c757d;">{{ project.desc }}</p>
              </div>
            </div>
          </section>

          <!-- 荣誉与奖项 -->
          <section>
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              荣誉与奖项
            </h2>
            <div class="ms-3">
              <div v-for="(award, idx) in data.awards" :key="idx" class="mb-3">
                <h3 class="h6 fw-semibold text-dark mb-0">{{ award.title }}</h3>
                <p class="text-muted mb-0 small">
                  {{ award.org }} · {{ award.year }}
                </p>
              </div>
            </div>
          </section>

          <!-- 联系我 -->
          <section>
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              联系我
            </h2>
            <div class="text-center">
              <div class="d-flex flex-wrap justify-content-center align-items-center gap-4 text-muted small">
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
          <p class="text-muted mb-0" style="font-size: 0.875rem; color: #6c757d;">
            <i class="fas fa-star me-1"></i>
            {{ data.footerText }}
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
      data: aboutData
    }
  }
}
</script>

<style scoped></style>

<!-- About.vue -->
<template>
  <div class="container py-4 px-3 view-container about-view">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8 typography-body">

        <div class="text-center mb-4">
          <h1 class="article-title mb-3">{{ pageTitle }}</h1>
        </div>

        <div class="d-flex flex-column" style="gap: 2rem;">

          <section>
            <h2 class="h4 fw-semibold pb-2 mb-3"
              style="border-bottom: 1px solid var(--app-custom-border-color); color: var(--app-text-emphasis);">
              {{ introductionTitle }}
            </h2>
            <p style="color: var(--app-text-secondary); margin-bottom: 0;">
              {{ data.introduction }}
            </p>
          </section>

          <section v-for="(entry, idx) in data.section" :key="idx">
            <h2 class="h4 fw-semibold pb-2 mb-3"
              style="border-bottom: 1px solid var(--app-custom-border-color); color: var(--app-text-emphasis);">
              {{ entry.title }}
            </h2>
            <div class="ms-3">
              <div v-for="(it, j) in entry.items" :key="j" class="mb-4">
                <p style="color: var(--app-text-muted); font-weight: bold; margin-bottom: 0.25rem;">{{ it.item }}</p>
                <p v-if="it.desc" style="color: var(--app-text-secondary); margin-bottom: 0;">{{ it.desc }}</p>
              </div>
            </div>
          </section>

          <section>
            <h2 class="h4 fw-semibold pb-2 mb-3"
              style="border-bottom: 1px solid var(--app-custom-border-color); color: 'var(--app-text-emphasis)';">
              {{ contactTitle }}
            </h2>
            <div class="text-center">
              <div
                style="display: flex; flex-wrap: wrap; justify-content: center; align-items: center; gap: 1rem; color: var(--app-text-secondary);">
                <a v-for="contact in data.contacts" :key="contact.label" :href="contact.link" target="_blank"
                  rel="noopener noreferrer" class="d-flex align-items-center text-decoration-none contact-link"
                  :style="{ color: 'var(--app-custom-icon-color)' }">
                  <i :class="contact.icon" class="me-1"></i>
                  <span>{{ contact.value }}</span>
                </a>
              </div>
            </div>
          </section>

        </div>

        <div class="text-center mt-5 pt-4" style="border-top: 1px solid var(--app-custom-border-color);">
          <p style="color: var(--app-text-secondary); margin-bottom: 0;">
            <i class="fas fa-star me-1"></i>
            {{ footerText }}
          </p>
        </div>

      </div>
    </div>
  </div>
</template>

<script>
import { useI18n } from 'vue-i18n'
import { loadAbout } from '@/utils/contentLoader'

/*
  AboutView
  - 关于页面
*/
export default {
  name: 'AboutView',
  setup() {
    const { t, locale } = useI18n()
    return { t, locale }
  },
  data() {
    return {
      data: {}
    }
  },
  computed: {
    pageTitle() {
      return this.t('about')
    },
    footerText() {
      return this.t('thanks')
    },
    introductionTitle() {
      return this.t('introduction')
    },
    contactTitle() {
      return this.t('contact')
    }
  },
  created() {
    this.loadAboutData()
  },
  watch: {
    locale() {
      this.loadAboutData()
    }
  },
  methods: {
    async loadAboutData() {
      try {
        this.data = await loadAbout() || {};
      } catch (error) {
        console.error('Failed to load about data:', error);
        this.data = {
          introduction: this.introductionTitle,
          section: [],
          contacts: []
        };
      }
    }
  }
}
</script>

<style scoped>
.typography-body {
  font-size: 1.125rem;
  line-height: 1.8;
}

.typography-body p {
  margin-bottom: 0.75rem;
}

.contact-link {
  transition: transform 0.2s ease;
}

.contact-link:hover {
  transform: translateY(-1px);
}
</style>
